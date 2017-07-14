"""
For looking at the impact of catheter order on the effective dose to
the DIL region in the prostate

INPUT: dwelltimes.txt, anatomy.txt, dvh_points.txt, report.txt
OUTPUT: BED to DIL when using optimized catheter order versus unoptimized
"""

import sys
import string
import itertools
import math
import numpy as np
import matplotlib.path as mpltPath
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay

if len(sys.argv) < 2: 
    print "Usage: python {} [input date]".format(sys.argv[0])

def read_dwelltimes(filename):
    '''Read dwelltimes.txt and return dwell positions and dwell times'''
    with open(filename) as f: 
        for i in range(5):
            # ignore unwanted lines
            f.readline()
        lines = f.readlines()

        coords = []
        times = []
        for line in lines[1:]:
            points = []
            x,y,z = map(float, line.split()[:3])
            points = [x,y,z]
            coords.append(points)
            time = map(float, line.split()[3:4])
            times.append(time)
        return coords, times

def read_points(filename):
    '''Reads dosepoint files, can be surface.txt, volumelow.txt, dvhpoints.txt etc.
    Returns the points for just the DIL region, can be changed to return just the points
    of the urethra as well'''
    organs = []
    with open(filename) as f:
        for i in range(4):
            f.readline()
            # ignore the "----" and "IPSA"
        coords = []
        organ_name = f.readline()
        num = f.readline()
        # collect the number listed after the organ name
        while True:
            line = f.readline()
            if line == "":
                # stop running if reached the end of file
                break
            if line[0] in string.digits + "-":
                # string.digits checks if string is a number
                coords.append(map(float, line.split()))
                # check if float and if so, append to coords
            else: # reached a new organ 
                # store the old one in the list and start a new organ 
                organs.append((organ_name, coords))
                coords = []
                organ_name = line.strip()
                num = f.readline()
    for name, coords in organs:
        if (name == "DIL"):
            # for this purpose only want the dose points that are in the DIL
            points = coords
    return points

def read_anatomy(filename):
    '''Read anatomy.txt and return the name of the organ and its coordinates'''
    organs = []
    with open(filename) as f:
        for i in range(4):
            f.readline()
            # ignore the "----" and "IPSA"
        coords = []
        organ_name = f.readline()
        num = f.readline()
        # collect the number listed after the organ name
        while True:
            line = f.readline()
            if line == "":
                # stop running if reached the end of file
                break
            if line[0] in string.digits + "-":
                # string.digits checks if string is a number
                coords.append(map(float, line.split()))
                # check if float and if so, add to coords
            else: # reached a new organ 
                # store the old one in the list and start a new organ 
                organs.append((organ_name, coords))
                coords = []
                organ_name = line.strip()
                num = f.readline()
    return organs

def read_report(filename):
    '''Reads report.txt and returns the volume of the DIL and the Target'''
    content = open(filename).read()
    DIL_vol = float(content.split("(DIL) volume = ")[1].split("\n")[0][:-1])
    # find where in the file it says "(DIL) volume = " and grab what comes after
    Target_vol = float(content.split("(Target) volume = ")[1].split("\n")[0][:-1])
    # repeat for the target
    return DIL_vol, Target_vol

def tetrahedron_volume(a, b, c, d):
    '''Used to determine whether point is inside organ; called by in_organ'''
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6

def convex_hull_volume(pts):
    '''Used to determine whether point is inside organ; called by in_organ'''
    ch = ConvexHull(pts)
    dt = Delaunay(pts[ch.vertices])
    tets = dt.points[dt.simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],tets[:, 2], tets[:, 3]))

def in_organ(coords, organs, name):
    '''Determines whether point is in the organ matching the name that is passed to the function
    Calls convex_hull_volume and tetrahedron_volume'''
    is_in = False
    for organ, o_coords in organs:
        if (name == organ):        
            organ_zs = []
            for line in o_coords:
                organ_z = line[2]
                organ_zs.append(organ_z)
            minz = min(organ_zs)
            maxz = max(organ_zs)
            if coords[2] < maxz and coords[2] > minz:
                organ_xs = []
                organ_ys = []
                for axes in o_coords:
                    organ_xs.append(axes[0])
                    organ_ys.append(axes[1])
                # turn the points into a numpy array instead of a list
                xpoints = np.array(organ_xs)
                ypoints = np.array(organ_ys)
                zpoints = np.array(organ_zs)
                # join all the coordinates into one array
                organ_points = np.vstack((xpoints, ypoints, zpoints)).T
                # find current volume of the points
                volume = convex_hull_volume(organ_points)
                
                # add the dwell point to the organ coordinates
                organ_xs.append(coords[0])
                organ_ys.append(coords[1])
                organ_zs.append(coords[2])
                # turn the new points into an array
                dwell_x = np.array(organ_xs)
                dwell_y = np.array(organ_ys)
                dwell_z = np.array(organ_zs)
                dwell_points = np.vstack((dwell_x, dwell_y, dwell_z)).T
                # find volume of new coordinates
                new_volume = convex_hull_volume(dwell_points)
                # if the two volumes are the same, that means the point is inside the organ
                if (new_volume == volume):
                    is_in = True
                else: 
                    is_in = False
            else:
                is_in = False
    return is_in

def convert_to_BED(doses, time, OPT=False):
    '''Converts the doses passed to BED. Takes the dose from all fractions as an argument but can be used for one fraction as well'''

    alpha = 0.15   #Gy^-1
    beta = 0.048   #Gy^-2
    mu = 61.6      #days^-1
    Tp = 42.0      #days
    Tk = 0.0       #days
    # convert all times from seconds to days
    duration = (time[1] / 0.75) / (24.0 * 3600.0) 
    D = 0.0
    N = len(doses)
    for i in range(N):
        # ensure that the doses passed to the function are in Gy
        D += doses[i]

    if (OPT):
        # if set to True, use optimized catheter placement
        # t1 is time spent in DIL
        t1 = (time[0]/ 0.75) / (24.0 * 3600.0)  
        # t2 is the difference between the full treatment time and the time spent in the DIL
        temp = time[1] - time[0]
        t2 = (temp/ 0.75) / (24.0 * 3600.0) 
        # the 0.75 factor accounts for additional time such as movement between catheters and lag time
        times = [t1, t2]
    else: 
        # default is to use one fraction of catheters is no particular order
        t = (time[1]/ 0.75) / (24.0 * 3600.0) 
        times = [t]

    g = 0.0
    for j in range(0,N):
        factor_1 = (doses[j]**2)*(times[j] - ((1/mu)*(1-math.exp(-mu*times[j]))))/(times[j]**2)
        factor_2 = 0.0
        for k in range(0,j):
            if (OPT):
                # time between fractions is the duration of the previous fraction
                factor_2 += (doses[k]*doses[j]*((math.exp(-mu*times[0]))*(math.exp(mu*times[k]) - 1)*(math.exp(-mu*times[j])-1)))/(times[k]*times[j])
            else:
                # there is no time between fractions when there is only one total fraction
                factor_2 += (doses[k]*doses[j]*((math.exp(-mu*0))*(math.exp(mu*times[k]) - 1)*(math.exp(-mu*times[j])-1)))/(times[k]*times[j])
        g += factor_1 - ((1/mu)*factor_2)
    g_factor = (2/mu)*g
    prolif = ((0.693/(alpha*Tp))*(duration-Tk))
    BED = (D + (1/((alpha/beta)) * g_factor)) - prolif
    # return dose in cGy
    return (BED*100.0)
           
def calculate_dose(sourcepoints, time, dosepoints):
    '''Calculates the dose using the TG-43 formulism and values for Ir-192 taken from IPSA'''
    
    D2 = ((dosepoints[0] - sourcepoints[0])*(dosepoints[0] -sourcepoints[0]) +
         (dosepoints[1] - sourcepoints[1])*(dosepoints[1] - sourcepoints[1]) +
         (dosepoints[2] - sourcepoints[2])*(dosepoints[2] - sourcepoints[2]))
    if (D2<=0.04):
        # cap distance at a minimum to avoid unreasonably large numbers
        D2=0.04
    D = math.sqrt(D2)

    # all values taken from structures.h in IPSA
    Activity = 8000.0
    AirKermaFactor = 4.03
    DoseRateConstant = 1.12
    AnisotropyFactor = 0.98
    HalfLife = 73.83
    g0 = 0.989054
    g1 = 0.0081319
    g2 = 0.0035177
    g3 = -0.00146637
    g4 = 0.000092437
    g5 = 0.0

    GeometryFactor = 1.0 / D2
    SourceStrength = AirKermaFactor * Activity
    RadialDoseFunction = g0 + g1*D + g2*D2 + g3*D2*D + g4*D2*D2 + g5*D2*D2*D

    doserate = SourceStrength * DoseRateConstant * AnisotropyFactor * RadialDoseFunction * GeometryFactor
    value = doserate / 3600.0
    Dose = value * time
    # return dose in cGy
    return Dose

def total_dose(dwell_info, organs, dvh_points):
    '''Main function, calculates the dose for each dose point and adds them together. Then converts the total physical dose for each point to BED'''

    DILtime, totaltime = total_time(dwell_info, organs)
    dwells = []
    # first loop to collect all points in the DIL
    for points, times in dwell_info:
        time = float(times[0])
        in_DIL = bool(in_organ(points, organs, "DIL"))
        dwells.append((in_DIL, points, time))
     
    numPoints = len(dvh_points)
    times = []
    times.append(DILtime)
    times.append(totaltime)
    total = 0.0
    complete_total = 0.0
    dvh_info = []
    for coords in dvh_points:
        # set all counters to 0
        inDIL = 0.0
        outDIL = 0.0
        complete_point = 0.0
        point_BED = 0.0
        for DIL, positions, time in dwells:
            dose = calculate_dose(positions, time, coords)
            # seperate the dose that comes from catheters inside and outside
            if (DIL):
                inDIL += dose
            else: 
                outDIL += dose
            # for the unoptimized version, collect all dose regardless of which catheter it comes from
            complete_point += dose
        # convert dose to Gy
        doses = [(inDIL/100.0), (outDIL/100.0)]
        complete_doses = [(complete_point/100.0)]
        # set OPT to True to optimize catheter positions
        point_BED = convert_to_BED(doses, times, OPT=True)
        # set OPT to false to calculate control
        complete_BED = convert_to_BED(complete_doses, times, OPT=False)
        
        # collect information that you want to print to .hist file
        dvh_info.append((coords, complete_BED))

        # total dose added across each point
        total += point_BED
        complete_total += complete_BED
    # average dose each dose point is getting
    average = total / numPoints
    average_complete = complete_total / numPoints 
    percent = ((total-complete_total)/complete_total) *100.0
    
    print "Average short: ", average
    print "Average long: ", average_complete
    print "Percent difference: ", percent 
    return dvh_info
    
def total_time(dwell_info, organs):
    ''' Calculates the total amount of time spent in the DIL region versus the total amount of time for the entire treatment'''
    total_time = 0.0
    DIL_time = 0.0
    for points, times in dwell_info:
        time = float(times[0])
        in_target = in_organ(points, organs, "Target")
        # if in the prostate, check if in DIL region
        if (in_target):
            # store the coordinates of the dwell position along with whether it is in
            # the DIL region and the dwell time
            in_DIL = bool(in_organ(points, organs, "DIL"))
            if (in_DIL):
                DIL_time += time
        total_time += time
    return DIL_time, total_time

def AddValueAtX(x, value, Histo, xmin, xmax, width):
    '''Add the value to the histogram'''
    numbin = 100.0
    # this code is mostly taken from IPSA and translated to python rather than c++
    if (x<xmin or x>xmax):
        return False
    xbin = math.floor((x-xmin)/width)
    if (xbin >= numbin):
        xbin = int(numbin - 1)
    else:
        xbin = int(xbin)
    Histo[xbin] += value
    return True

def fill_DVH(dvh_info):
    '''Create the histogram file for either the optimized or unoptimized treatment'''
    # this code is mostly taken from IPSA and translated to python rather than c++
    xmax = 0.0
    doses = []
    for coords, dose in dvh_info:
        doses.append(dose)
        if (dose > xmax):
                xmax = dose
    xmin = 0.0
    width = (xmax - xmin)/100.0
    Histo = [0] * 100
    gridSpacing = 0.2
    sliceZ = []
    numPoints = len(dvh_points)

    xs = []
    ys = []
    zs = []    
    for coords in dvh_points:
        xs.append(coords[0])
        ys.append(coords[1])
        zs.append(coords[2])

    for i in range(0,numPoints):
        if not sliceZ:
            sliceZ.append(zs[i])
        elif zs[i] != sliceZ[-1]:
            sliceZ.append(zs[i])

    sliceThickness = []
    if (len(sliceZ) < 2):
        print "Warning, less than 2 slices. Quitting!"
        return
    else:
        # first and last slices are treated outside of the loop
        sliceThickness.append(abs(sliceZ[0] - sliceZ[1]))
        for j in range(1, len(sliceZ)-1):
            zDiff1 = abs(sliceZ[j] - sliceZ[j-1])
            zDiff2 = abs(sliceZ[j] - sliceZ[j+1])
            sliceThickness.append(0.5*zDiff1 + 0.5*zDiff2) 
        sliceThickness.append(abs(sliceZ[-1] - sliceZ[-2]))
                              
    if (len(sliceZ) != len(sliceThickness)):
        print "sliceZ: ", len(sliceZ)
        print "sliceThickness: ", len(sliceThickness)
        print "Warning, sliceZ and sliceThickness are different sizes. Quitting!"
        return
                    
    z_previousPoint = 999999999999999999.9
    organVolume = 0.0
    sliceZIndex = 9999999999
    sliceVolume = [0] * len(sliceZ)
                        
    for iPoint in range(0, numPoints):
        getNewThicknessIndex = False
        z_thisPoint = zs[iPoint]
        if (iPoint == 0):
            getNewThicknessIndex = True
        elif (z_thisPoint != z_previousPoint):
            getNewThicknessIndex = True
        if (getNewThicknessIndex):
            for iSlice in range(0,len(sliceZ)):
                if (z_thisPoint == sliceZ[iSlice]):
                    sliceZIndex = iSlice
                              
        volumeElement = sliceThickness[sliceZIndex] * gridSpacing * gridSpacing
        sliceVolume[sliceZIndex] += volumeElement
        organVolume += volumeElement
        measuredDose = doses[iPoint]
        AddValueAtX(measuredDose, volumeElement, Histo, xmin, xmax, width)
        z_previousPoint = z_thisPoint
                              
    Data = []
    xbin  = 0
    for number in Histo:
        x = xmin + (xbin * width)
        value = number
        Data.append((x, value))
        xbin += 1
                              
    newfile = open("DIL_DVH.hist", "w")
    newfile.write("\n".join(" ".join(map(str, r)) for r in Data))
                            
# Read in the values from the input files
coords, times = read_dwelltimes(sys.argv[1])
organs = read_anatomy(sys.argv[2])
dvh_points = read_points(sys.argv[3])
DILvol, Tarvol = read_report(sys.argv[4])

# Zip the coordinates and the times for the dwells together 
dwell_info = zip(coords, times)

# Call the main function and fill the DVH
dvh_info = total_dose(dwell_info, organs, dvh_points)   
fill_DVH(dvh_info)


