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

staveTop.ScaleImage(1)
staveBottom.ScaleImage(1)

#get the scaled image
img_edges = staveTop.getImage()

#finding the staves - triggers an algorithm that looks for the stave, using relative coordinates
staveTop.FindStaveWithin(0,1.0,0,0.5)
staveBottom.FindStaveWithin(0,1.0,0.5,1.0)

#print the positions of the staves
staveTop.Echo()
staveBottom.Echo()

#create a deep copy of the image, to which the edges/regions will be drawn
staveTop.DrawEdges(img_edges)
staveBottom.DrawEdges(img_edges)


numModules = 14

#large regions
for i in range(numModules):
  staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.0,0.5,"large")
  staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.0,0.5,"large")

#large regions - return pipe
for i in reversed(range(numModules)):
  staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.5,1.0,"large")
  staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.5,1.0,"large")

#small regions above the pipe
for i in range(numModules):
  staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.28260-0.03478,0.28260+0.03478,"small")
  staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.28260-0.03478,0.28260+0.03478,"small")

#small regions above the pipe (return pipe)
for i in reversed(range(numModules)):
  staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.5+0.28260-0.03478,0.5+0.28260+0.03478,"small")
  staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.5+0.28260-0.03478,0.5+0.28260+0.03478,"small")

#drawing the regions
staveBottom.DrawRegions(img_edges,"large")
staveTop.DrawRegions(img_edges,"large")
staveBottom.DrawRegions(img_edges,"small")
staveTop.DrawRegions(img_edges,"small")

staveBottomTemp = staveBottom.getTemperatures("small")
staveTopTemp = staveTop.getTemperatures("small")



#using the temperature profile provided by the FEA (otherwise a linear extrapolation would be used)
temperatureProfile = [0,0.062,0.114,0.152,0.188,0.224,0.26,0.295,0.33,0.364,0.398,0.432,0.466,0.499,0.533,0.568,0.603,0.637,0.672,0.706,0.74,0.774,0.807,0.841,0.873,0.906,0.937,0.969,1.0]
staveBottom.setTemperatureProfile(temperatureProfile)
staveTop.setTemperatureProfile(temperatureProfile)

#extracting the impedances
largeTop = staveTop.getImpedances("large")
largeBottom = staveBottom.getImpedances("large")
smallTop = staveTop.getImpedances("small")
smallBottom = staveBottom.getImpedances("small")

#adding the small region impedances from top and bottom in parallel:
smallCombined = list(1/(1/np.array(smallTop) + 1/np.array(smallBottom)))

#plotting
plt.figure(figsize=(12,6))
plt.plot(largeTop, label="Large Region: top")
plt.plot(largeBottom, label="Large Region: bottom")
plt.plot(smallTop, label="Small Region: top")
plt.plot(smallBottom, label="Small Region: bottom")
plt.plot(smallCombined, label="Small Region: combined")
plt.xlabel("Module number")
plt.ylabel("Thermal Impedance [C/W]")
plt.title("Thermal Impedances for the Stave control regions")
yrange = int(1+1.1*np.max([np.max(largeTop),np.max(largeBottom),np.max(smallTop),np.max(smallBottom)]))
plt.xticks(np.arange(0, 28, 1.0))
plt.yticks(np.arange(0, yrange, 0.5))
plt.axis([-0.5,27.5,0,yrange])
plt.grid()
plt.legend()
plt.show()

"""
plt.imshow(img_edges)
plt.show()
"""