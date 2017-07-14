"""Script to run draw2Histograms.py over all the organs in the files directory of 2 cases"""

import sys
import os
import draw2Histograms

dir_1 = sys.argv[1]
dir_2 = sys.argv[2]

for file in os.listdir(dir_1):
    print dir_1, dir_2, file
    # draw histograms for any file that ends in .hist and has a partner in the other directory
    if file.endswith(".hist"):
        drawHistogram2.drawHistogram(dir_1 + "/" + file, dir_2 + "/" + file)

