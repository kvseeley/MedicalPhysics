"""Draws a histogram given a .hist file"""

import os
import sys
from ROOT import TH1F, TCanvas

canvas = TCanvas()

def drawHistogram(file_name):
    print "Reading from " + file_name 

    xs = [ ]
    ys = [ ]

    with open(file_name) as f: 
        for line in f.readlines():
            x,y = map(float, line.split("  "))
            xs.append(x)
            ys.append(y)

    xmin = min(xs)
    xmax = max(xs)

    points = zip(xs,ys)

    base_name = os.path.basename(file_name)
    no_extension = base_name.split(".")[0] 

    if no_extension.endswith("Cum"):
        name = no_extension[:-3] + " Cumulative DVH"
    else:
        name = no_extension + " DVH"

    hist = TH1F(name, name, len(points), xmin, xmax)
   
    for x,y in points:
       hist.Fill(x,y)

    hist.GetXaxis().SetTitle("Dose (cGy)")
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().SetTitle("Volume (cc)")
    hist.GetYaxis().CenterTitle()

    #canvas = TCanvas()
    hist.Draw()

    #canvas.Update()
    canvas.SaveAs(no_extension+"DVH.png")

    #return canvas

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "drawHistogram.py needs a file name as an argument"
        exit(1)
    drawHistogram(sys.argv[1])
    raw_input ("Press Enter to close")
