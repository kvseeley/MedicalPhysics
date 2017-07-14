"""Given an HDR Prostate case, can change the dwells files in order to resemble a PPI case"""

import sys
import itertools
import math
import numpy as np

if len(sys.argv) < 2: 
    print "Usage: python {} [input date]".format(sys.argv[0])

def read(filename):
    '''Read the file and return the max and min x,y,z values'''
    with open(filename) as f:
        lines = f.readlines()

    zcoord = 0
    circles = []
    pairs = []
    coords = []
    zs = []
    i = 0
    for line in lines[1:]:
        x,y,z = map(float, line.split())
        coords = [x,y]
        zs.append(z)
        # if first iteration then just add the z coord
        if i == 0:
            pairs.append(coords)
            zcoord = z
            i = 1
        elif z == zcoord:
            pairs.append(coords)
        # if z coords do not match, have not moved on to a new circle
        else: 
            circles.append(pairs)
            pairs = []
            zcoord = z

    circles.remove(circles[0])
    x = []
    y = []

    for pairs in circles: 
        for coords in pairs:
            x.append(coords[0])
            y.append(coords[1])

    xmax = max(x)
    xmin = min(x)
    ymax = max(y)
    ymin = min(y)
    zmax = max(zs)
    zmin = min(zs)
    return xmin, xmax, ymin, ymax, zmin, zmax

def grid(xmin, xmax, ymin, ymax, zmin, zmax):
    '''Creates the grid and prints it to a file'''

    # define starting positions
    xstart = math.floor(xmin)
    xstop = math.ceil(xmax)
    ystart = math.floor(ymin)
    ystop = math.ceil(ymax)
    zstart = math.floor(zmin)
    zstop = math.ceil(zmax)
    
    xspace = np.abs(xstop - xstart) / 0.5
    yspace = np.abs(ystop - ystart) / 0.5
    zspace = np.abs(zstop - zstart) / 0.5
    xs = np.linspace(xstart, xstop, xspace, endpoint=False)
    ys = np.linspace(ystart, ystop, yspace, endpoint=False)
    zs = np.linspace(zstart, zstop, zspace, endpoint=False)
   
    positions = []
    rows = []
    oldy = 0
    i = 0
    needle = 0
    for x in xs:
        for y in ys:
            for z in zs:
                positions.append(x)
                positions.append(y)
                positions.append(z)
                if i == 0:
                    oldy = y
                    # additional values needed for the dwells file
                    positions.extend([needle, "inactive", "unfrozen"])
                    rows.append(positions)
                    positions = []
                    i = 1
                else:
                    # if the y values are not the same, has moved to another needle
                    if y != oldy:
                        needle += 1
                        oldy = y
                    positions.extend([needle, "inactive", "unfrozen"])
                    rows.append(positions)
                    positions = []
    
    # write the file
    newfile = open("new_dwells.txt", "w")
    newfile.write("_______________________________________________________________________\n")
    newfile.write("                                  IPSA                                 \n")
    newfile.write("_______________________________________________________________________\n")
    newfile.write("\n".join(" ".join(map(str, r)) for r in rows))
    return rows

# Read in the values and run the main function
xmin, xmax, ymin, ymax, zmin, zmax = read(sys.argv[1])
newgrid = grid(xmin, xmax, ymin, ymax, zmin, zmax)


