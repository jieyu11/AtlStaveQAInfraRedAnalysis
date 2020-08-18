#!/usr/bin/env python

'''
impedanceFromCSV.py

Author: Lubos Vozdecky (Queen Mary, University of London)
About: This program takes a thermal image of an ATLAS Itk Stave Support in CSV format,
  finds the stave and computes the thermal impedances using data of the cooling fluid saved in parameters.cfg.
  The program outputs the result impedances into the /data folder as a CSV file.
'''

import numpy as np
import os
import csv
import logging
import argparse
from matplotlib import pyplot as plt
from stave import Stave


parser = argparse.ArgumentParser()
parser.add_argument("path", help="The path to the input CSV file")
parser.add_argument("config", help="The path to the configuration file")
parser.add_argument("-d","--debug", help="Runs the code in debug mode", action="store_true")
parser.add_argument("-g","--graphs", help="Outputs the graph", action="store_true")
parser.add_argument("-1f","--one_face", help="Using IR image with one face only", action="store_true")
args = parser.parse_args()


inputFile = args.path
configFile = args.config

#check if the suffix is .csv
if inputFile[-3:] != "csv":
  print("The input file should be a .csv file")
  quit()

#delete the debug folder if it exists and create a new one
if args.debug:
  if "debug_output" in os.listdir("."):
    os.system("rm -r debug_output")
  os.mkdir("debug_output")

#create the output folder if it doesn't exist
if not "output" in os.listdir("."):
  os.system("mkdir output")

#set up the debugging log
if args.debug:
  logging.basicConfig(filename='debug_output/debug.log',level=logging.DEBUG)


gitHash = os.popen('git rev-parse --short HEAD').read()[:-2]
gitDate = os.popen('git log -1 --format=%cd').read()

logging.debug("Running the code version: " + gitHash + " " + gitDate)

#fetch the CSV file
logging.debug("Opening the CSV file")
imgList = []
with open(inputFile) as csvfile:
  reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    imgList.append(row)
image = np.array(imgList)

#creating the staves + loading the parameters from the config file
staveTop = Stave(image, configFile)
if not(args.one_face):
  staveBottom = Stave(image, configFile)

#scale up the images
staveTop.ScaleImage(10)
if not(args.one_face):
  staveBottom.ScaleImage(10)

#get the scaled image
img_edges = staveTop.getImage()

#finding the staves - triggers an algorithm that looks for the stave, using relative coordinates
if args.one_face:
  staveTop.FindStaveWithin(0,1.0,0,1.0)
else:
  staveTop.FindStaveWithin(0,1.0,0,0.5)
  staveBottom.FindStaveWithin(0,1.0,0.5,1.0)
  
#print the positions of the staves
print("Staves' edges found at:")
staveTop.Echo()
if not(args.one_face):
  staveBottom.Echo()

#create a deep copy of the image, to which the edges/regions will be drawn
staveTop.DrawEdges(img_edges)
if not(args.one_face):
  staveBottom.DrawEdges(img_edges)


numModules = 14

#large regions
for i in range(numModules):
  staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.0,0.5,"large")
  if not(args.one_face):
    staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.0,0.5,"large")

#large regions - return pipe
for i in reversed(range(numModules)):
  staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.5,1.0,"large")
  if not(args.one_face):
    staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.5,1.0,"large")

#small regions above the pipe
for i in range(numModules):
  #exception for near-edge regions
  if i == 0:
    staveTop.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
    if not(args.one_face):
      staveBottom.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
  elif i==13:
    staveTop.AddUBendRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.0174545,0.247826,0.317391,0.13,0.0869565,"small",bend="downwards")
    if not(args.one_face):
      staveBottom.AddUBendRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.0174545,0.247826,0.317391,0.13,0.0869565,"small",bend="downwards")
  else:
    staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
    if not(args.one_face):
      staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
    
#small regions above the pipe (return pipe)
for i in reversed(range(numModules)):
  if i == 0:
    staveTop.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
    if not(args.one_face):
      staveBottom.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
  elif i==13:
    staveTop.AddUBendRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.0174545,0.682609,0.752174,0.13,0.0869565,"small",bend="upwards")
    if not(args.one_face):
      staveBottom.AddUBendRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.0174545,0.682609,0.752174,0.13,0.0869565,"small",bend="upwards")
  else:
    staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
    if not(args.one_face):
      staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")

#drawing the regions

staveTop.DrawRegions(img_edges,"large")
staveTop.DrawRegions(img_edges,"small")

staveTopTemp = staveTop.getTemperatures("small")


if not(args.one_face):
  staveBottom.DrawRegions(img_edges,"small")
  staveBottom.DrawRegions(img_edges,"large")

  staveBottomTemp = staveBottom.getTemperatures("small")


#extracting the impedances
largeTop = staveTop.getImpedances("large")
smallTop = staveTop.getImpedances("small")

if not(args.one_face):
  largeBottom = staveBottom.getImpedances("large")
  smallBottom = staveBottom.getImpedances("small")

#compute the combined impedance
smallTopThere = np.array(smallTop[0:14])
smallTopReturn = np.array(smallTop[14:28])
impedanceCombinedTop = 1/(1/smallTopThere + 1/np.flip(smallTopReturn))

if not(args.one_face):
  smallBottomThere = np.array(smallBottom[0:14])
  smallBottomReturn = np.array(smallBottom[14:28])
  impedanceCombinedBottom = 1/(1/smallBottomThere + 1/np.flip(smallBottomReturn))

#savign data into the CSV file
outputFilename = "output/" + inputFile.split("/")[-1][:-4] + "_IMPEDANCES"
print("Outputing data into a file: " + outputFilename + ".csv")
with open(outputFilename+".csv", "w+") as f:
  if args.one_face:
    f.write('#, topLargeRegion, topSmallRegion, smallRegionCombinedTop \n')
  else:
    f.write("#, topLargeRegion, bottomLargeRegion, topSmallRegion, bottomSmallRegion, smallRegionCombinedTop, smallRegionCombinedBottom \n")
  for i in range(0,28):
    if i<14:
      if args.one_face:
        f.write(str(i)+", "+str(largeTop[i])+", "+str(smallTop[i])+", "+str(impedanceCombinedTop[i]) + "\n")
      else:
        f.write(str(i)+", "+str(largeTop[i])+", "+str(largeBottom[i])+", "+str(smallTop[i])+", "+str(smallBottom[i]) + ", "+str(impedanceCombinedTop[i]) + ", "+str(impedanceCombinedBottom[i]) + "\n")
    else:
      if args.one_face:
        f.write(str(i)+", "+str(largeTop[i])+", "+str(smallTop[i])+"\n")
      else:
        f.write(str(i)+", "+str(largeTop[i])+", "+str(largeBottom[i])+", "+str(smallTop[i])+", "+str(smallBottom[i]) + "\n")
f.close()


#plotting if -g option selected
if args.graphs:
  plt.figure(figsize=(12,6))
  plt.plot(largeTop, label="Large Region: top")
  plt.plot(smallTop, label="Small Region: top")
  plt.plot(impedanceCombinedTop, label="Small Region: top combined")
  if not(args.one_face):
    plt.plot(largeBottom, label="Large Region: bottom")
    plt.plot(smallBottom, label="Small Region: bottom")
    plt.plot(impedanceCombinedBottom, label="Small Region: bottom combined")
  plt.xlabel("Region number")
  plt.ylabel("Thermal Impedance [K/W]")
  plt.title("Thermal Impedances for the Stave control regions")
  if args.one_face:
    yrange = int(1+1.1*np.max([np.max(largeTop),np.max(smallTop)]))
  else:
    yrange = int(1+1.1*np.max([np.max(largeTop),np.max(largeBottom),np.max(smallTop),np.max(smallBottom)]))
  plt.xticks(np.arange(0, 28, 1.0))
  plt.yticks(np.arange(0, yrange, 0.5))
  plt.axis([-0.5,27.5,0,yrange])
  plt.grid()
  plt.legend()
  plt.text(0, -0.13*yrange, "Code version: " + gitHash + " " + gitDate[:-6], fontsize=10)
  plt.savefig(outputFilename + ".png")
  print("Outputing graphical output into a file: " + outputFilename + ".png")
  

if args.debug:
  plt.clf()
  plt.imshow(img_edges)
  plt.savefig("debug_output/edges.png")
