#!/usr/bin/env python
import numpy as np
import cv2
#import ROOT
import csv
import logging
import argparse
from matplotlib import pyplot as plt
from stave import Stave


parser = argparse.ArgumentParser()
parser.add_argument("path", help="The path to the input CSV file")
parser.add_argument("-d","--debug", help="Runs the code in debug mode", action="store_true")
parser.add_argument("-g","--graphs", help="Outputs the graph", action="store_true")
args = parser.parse_args()

inputFile = args.path

#check if the suffix is .csv
if inputFile[-3:] != "csv":
  print("The input file should be a .csv file")
  quit()

#create the folder for debug output
try:
  if args.debug: 
	os.mkdir("debug_output")
	print("Running the code in the ***DEBUG MODE***")
except OSError:
  print("Folder called 'debug_output' already exists. Please delete it or rename it.")
  quit()


#fetch the CSV file
logging.debug("Opening the CSV file")
imgList = []
with open(inputFile) as csvfile:
  reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    imgList.append(row)
image = np.array(imgList)

#creating the staves + loading the parameters from the config file
staveTop = Stave(image, "parameters.cfg",args)
staveBottom = Stave(image, "parameters.cfg",args)

height = image.shape[0]
width = image.shape[1]
#finding the staves - triggers an algorithm that looks for the stave
staveTop.FindStaveWithin(0,width,0,height/2)
staveBottom.FindStaveWithin(0,width,height/2,height)

staveTop.Echo()
staveBottom.Echo()

img_edges = np.copy(image)
staveTop.DrawEdges(img_edges)
staveBottom.DrawEdges(img_edges)

plt.imshow(img_edges)
plt.show()


