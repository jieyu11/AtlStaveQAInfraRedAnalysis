#!/usr/bin/env python

'''
impedanceFromCSV.py

Author: Lubos Vozdecky (Queen Mary, University of London)
About: This program takes a thermal image of an ATLAS Itk Stave Support in CSV format,
  finds the stave and computes the thermal impedances using data of the cooling fluid saved in parameters.cfg.
  The program outputs the result impedances into the /data folder as a CSV file.
'''

import numpy as np
import cv2
import os
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

#set up the debugging log
if args.debug:
  logging.basicConfig(filename='debug_output/debug.log',level=logging.DEBUG)

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

#scale up the images
staveTop.ScaleImage(10)
staveBottom.ScaleImage(10)

#get the scaled image
img_edges = staveTop.getImage()

#finding the staves - triggers an algorithm that looks for the stave, using relative coordinates
staveTop.FindStaveWithin(0,1.0,0,0.5)
staveBottom.FindStaveWithin(0,1.0,0.5,1.0)

#print the positions of the staves
print("Staves' edges found at:")
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
  #exception for near-edge regions
  if i == 0:
    staveTop.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
    staveBottom.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
  elif i==13:
    staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.5/numModules,0.247826,0.317391,"small")
    staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.5/numModules,0.247826,0.317391,"small")
  else:
    staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
    staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
    
#small regions above the pipe (return pipe)
for i in reversed(range(numModules)):
  if i == 0:
    staveTop.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
    staveBottom.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
  elif i==13:
    staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.5/numModules,0.682609,0.752174,"small")
    staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.5/numModules,0.682609,0.752174,"small")
  else:
    staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
    staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")

#drawing the regions
staveBottom.DrawRegions(img_edges,"large")
staveTop.DrawRegions(img_edges,"large")
staveBottom.DrawRegions(img_edges,"small")
staveTop.DrawRegions(img_edges,"small")

staveBottomTemp = staveBottom.getTemperatures("small")
staveTopTemp = staveTop.getTemperatures("small")



#using the temperature profile provided by the FEA (otherwise a linear extrapolation would be used)
#temperatureProfile = [0,0.062,0.114,0.152,0.188,0.224,0.26,0.295,0.33,0.364,0.398,0.432,0.466,0.499,0.533,0.568,0.603,0.637,0.672,0.706,0.74,0.774,0.807,0.841,0.873,0.906,0.937,0.969,1.0]
#FEA:
temperatureProfile = [0.000,0.060,0.109,0.143,0.178,0.213,0.247,0.281,0.316,0.350,0.384,0.418,0.451,0.485,0.520,0.555,0.589,0.623,0.657,0.691,0.725,0.760,0.794,0.828,0.862,0.896,0.930,0.965,1.000]

staveBottom.setTemperatureProfile(temperatureProfile)
staveTop.setTemperatureProfile(temperatureProfile)

#extracting the impedances
largeTop = staveTop.getImpedances("large")
largeBottom = staveBottom.getImpedances("large")
smallTop = staveTop.getImpedances("small")
smallBottom = staveBottom.getImpedances("small")

#compute the combined impedance
smallTopThere = np.array(smallTop[0:14])
smallTopReturn = np.array(smallTop[14:28])
impedanceCombinedTop = 1/(1/smallTopThere + 1/smallTopReturn)

smallBottomThere = np.array(smallBottom[0:14])
smallBottomReturn = np.array(smallBottom[14:28])
impedanceCombinedBottom = 1/(1/smallBottomThere + 1/smallBottomReturn)

#savign data into the CSV file
outputFilename = "output/" + inputFile.split("/")[-1][:-4] + "_IMPEDANCES"
print("Outputing data into a file: " + outputFilename + ".csv")
with open(outputFilename+".csv", "w+") as f:
  f.write("#, topLargeRegion, bottomLargeRegion, topSmallRegion, bottomSmallRegion, smallRegionCombinedTop, smallRegionCombinedBottom \n")
  for i in range(0,28):
    if i<14:
      f.write(str(i)+", "+str(largeTop[i])+", "+str(largeBottom[i])+", "+str(smallTop[i])+", "+str(smallBottom[i]) + ", "+str(impedanceCombinedTop[i]) + ", "+str(impedanceCombinedBottom[i]) + "\n")
    else:
      f.write(str(i)+", "+str(largeTop[i])+", "+str(largeBottom[i])+", "+str(smallTop[i])+", "+str(smallBottom[i]) + "\n")
f.close()


#plotting if -g option selected
if args.graphs:
  plt.figure(figsize=(12,6))
  plt.plot(largeTop, label="Large Region: top")
  plt.plot(largeBottom, label="Large Region: bottom")
  plt.plot(smallTop, label="Small Region: top")
  plt.plot(smallBottom, label="Small Region: bottom")
  plt.plot(impedanceCombinedTop, label="Small Region: top combined")
  plt.plot(impedanceCombinedBottom, label="Small Region: bottom combined")
  plt.xlabel("Region number")
  plt.ylabel("Thermal Impedance [K/W]")
  plt.title("Thermal Impedances for the Stave control regions")
  yrange = int(1+1.1*np.max([np.max(largeTop),np.max(largeBottom),np.max(smallTop),np.max(smallBottom)]))
  plt.xticks(np.arange(0, 28, 1.0))
  plt.yticks(np.arange(0, yrange, 0.5))
  plt.axis([-0.5,27.5,0,yrange])
  plt.grid()
  plt.legend()
  plt.savefig(outputFilename + ".png")
  print("Outputing graphical output into a file: " + outputFilename + ".png")
  

if args.debug:
  plt.clf()
  plt.imshow(img_edges)
  plt.savefig("debug_output/edges.png")
