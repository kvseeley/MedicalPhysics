"""Draws histogram for several organs in a directory"""

import os
from drawHistogram import drawHistogram 
from ROOT import TCanvas, TPad

canvas = TCanvas()
canvas.Divide(2,2)

canvas.cd(0)
drawHistogram("outputs/output/files/Rectum.hist")
canvas.cd(1)
drawHistogram("outputs/output/files/Bladder.hist")
canvas.cd(2)
drawHistogram("outputs/output/files/Urethra.hist")
canvas.cd(3)
drawHistogram("outputs/output/files/Target.hist")


