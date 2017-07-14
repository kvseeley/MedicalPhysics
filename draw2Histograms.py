"""Draws two histograms and overlays them on the same canvas given two .hist files"""

import os
import sys
from ROOT import TH1F, TCanvas, THistPainter, TLegend

canvas = TCanvas()
def drawHistogram(file_name, file_name2):
    print "Reading from " + file_name + " and " + file_name2
    xs = [ ]
    ys = [ ]
    xs2 = [ ]
    ys2 = [ ]

    # fill histogram one coordinates
    with open(file_name) as f: 
        for line in f.readlines():
            x,y = map(float, line.split())
            xs.append(x)
            ys.append(y)
    xmin = min(xs)
    xmax = max(xs)
    points = zip(xs,ys)

    # fill histogram two coordinates
    with open(file_name2) as f2: 
        for line in f2.readlines():
            x2,y2 = map(float, line.split())
            xs2.append(x2)
            ys2.append(y2)
    xmin2 = min(xs2)
    xmax2 = max(xs2)
    points2 = zip(xs2,ys2)

    # use the information in the filename to create the name of the histograms
    name = " "
    base_name = os.path.basename(file_name)
    no_extension = base_name.split(".")[0] 

    if no_extension.endswith("Cum"):
        name = no_extension[:-3] + " Cumulative DVH"
    else:
        name = no_extension + " Natural DVH"
    base_name2 = os.path.basename(file_name2)
    no_extension2 = base_name2.split(".")[0] 

    if no_extension2.endswith("Cum"):
        name2 = no_extension2[:-3] + " Cumulative DVH"
    else:
        name2 = no_extension2 + " DVH"
    maxtotal = max(xs + xs2)
    mintotal = min(xs + xs2)

    # create the two histograms
    hist = TH1F(name, name, len(points), mintotal, maxtotal)
    hist2 = TH1F(name2, name2, len(points), mintotal, maxtotal)

    # fill the two histograms
    for x,y in points:
        hist.Fill(x,y)
    for x2,y2 in points2:
        hist2.Fill(x2,y2)

    # set title and style for histograms
    hist.GetXaxis().SetTitle("Dose (cGy)")
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().SetTitle("Volume (cc)")
    hist.GetYaxis().CenterTitle()
    # change line color for the second histogram (Red)
    hist2.SetLineColor(2)

    # draw the histograms on the same canvas
    hist.Draw()
    hist2.Draw("same")

    # create a legend with the coordinates of where you want it located
    legend = TLegend(0.75, 0.5, 0.97, 0.7)
    legend.SetTextSize(0.035)
    # legend names need to be changed for each study
    legend.AddEntry(hist, "Short Treatment")
    legend.AddEntry(hist2, "Long Treatment")
    legend.Draw()

    # save the histogram as a .png file
    canvas.SaveAs(no_extension+".png")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "drawHistogram.py needs a file name as an argument"
        exit(1)
    drawHistogram(sys.argv[1],sys.argv[2])
    raw_input ("Press Enter to close")
