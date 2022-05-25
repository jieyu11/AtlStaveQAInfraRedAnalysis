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
import configparser
from matplotlib import pyplot as plt
from stave import Stave

parser = argparse.ArgumentParser()
parser.add_argument("path", help="The path to the input CSV file")
parser.add_argument("config", help="The path to the configuration file")
parser.add_argument("-d","--debug", help="Runs the code in debug mode", action="store_true")
parser.add_argument("-g","--graphs", help="Outputs the graph", action="store_true")
parser.add_argument("-1f","--one_face", help="Using IR image with one face only", action="store_true")
parser.add_argument("-L","--L_flip", help="Flips the input image horizontally  before processing.", action="store_true")
parser.add_argument("-J","--J_flip", help="Flips the input image both vertically and horizontally before processing.", action="store_true")
args = parser.parse_args()


inputFile = args.path
configFile = args.config

#check if the suffix is .csv
if inputFile[-3:] != "csv":
  print("The input file should be a .csv file")
  quit()

#delete the debug folder if it exists and create a new one
if args.debug:
  if not os.path.isdir("debug_output"):
    os.mkdir("debug_output")
  if os.path.isfile("debug_output/debug.log"):
    os.remove("debug_output/debug.log")

#create the output folder if it doesn't exist
if not "output" in os.listdir("."):
  os.system("mkdir output")

#set up the debugging log
if args.debug:
  logging.basicConfig(filename='debug_output/debug.log',level=logging.DEBUG)

if args.L_flip and args.J_flip:
  print("Cannot have both L and J flip options activated at the same time. Exiting...")
  exit()

#get the git version of the code so it can be printed on the output graphs
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

if args.L_flip:
  image = np.flip(image,axis=1)

if args.J_flip:
  image = np.flip(image,axis=0)
  image = np.flip(image,axis=1)

#creating the staves + loading the parameters from the config file
staveTop = Stave(image, configFile)
if not(args.one_face):
  staveBottom = Stave(image, configFile)

#scale up the images with linear extrapolation to get better results for small regions
staveTop.ScaleImage(10)
if not(args.one_face):
  staveBottom.ScaleImage(10)

#get the scaled image; it's used later for showing the regions for debugging
img_edges = staveTop.getImage()

#finding the staves - triggers an algorithm that looks for the stave, using relative coordinates
if args.one_face:
  staveTop.FindStaveWithin(0,1.0,0,1.0)
else:
  staveTop.FindStaveWithin(0,1.0,0,0.48)
  staveBottom.FindStaveWithin(0,1.0,0.54,1.0)
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
    staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.5,1.0,"large")

#large regions - return pipe
for i in reversed(range(numModules)):
  staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.5,1.0,"large")
  if not(args.one_face):
    staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.0,0.5,"large")

#small regions above the pipe
for i in range(numModules):
  #exception for near-edge regions
  if i == 0:
    staveTop.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
    if not(args.one_face):
      staveBottom.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
  elif i==13:
    staveTop.AddUBendRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.0174545,0.247826,0.317391,0.13,0.0869565,"small",bend="downwards")
    if not(args.one_face):
      staveBottom.AddUBendRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.0174545,0.682609,0.752174,0.13,0.0869565,"small",bend="upwards")
  else:
    staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
    if not(args.one_face):
      staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
#small regions above the pipe (return pipe)
for i in reversed(range(numModules)):
  #exception for near-edge regions
  if i == 0:
    staveTop.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
    if not(args.one_face):
      staveBottom.AddRegion(i*1.0/numModules + 0.1/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")
  elif i==13:
    staveTop.AddUBendRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.0174545,0.682609,0.752174,0.13,0.0869565,"small",bend="upwards")
    if not(args.one_face):
      staveBottom.AddUBendRegion(i*1.0/numModules,(i+1)*1.0/numModules - 0.0174545,0.247826,0.317391,0.13,0.0869565,"small",bend="downwards")
  else:
    staveTop.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.682609,0.752174,"small")
    if not(args.one_face):
      staveBottom.AddRegion(i*1.0/numModules,(i+1)*1.0/numModules,0.247826,0.317391,"small")

#end-of-stave ear region
#the region is defined to stay safely away from the edges: in x direction 0.1 of module length is subtracted from both sides; for y direction it's 5% of the stave width
staveTop.AddRegion(0.1/numModules,154.0/1375-0.1/numModules,-49.0/115+0.05,0.0,"ear")
if not args.one_face:
  staveBottom.AddRegion(0.1/numModules,154.0/1375-0.1/numModules,1.0,1.0+49.0/115-0.05,"ear")

#drawing the regions
staveTop.DrawRegions(img_edges,"large")
staveTop.DrawRegions(img_edges,"small")
staveTop.DrawRegions(img_edges,"ear")

staveTopTemp = staveTop.getTemperatures("small")


if not(args.one_face):
  staveBottom.DrawRegions(img_edges,"small")
  staveBottom.DrawRegions(img_edges,"large")
  staveBottom.DrawRegions(img_edges,"ear")

  staveBottomTemp = staveBottom.getTemperatures("small")

#correcting the temperature for the regions around the EoS ear (see Documents/2020-09-09-EOS-Impedances.pdf)
logging.debug("Importing variables from config file " + configFile + " in the impedanceFromCSV.py script.")
config = configparser.ConfigParser()
config.read(configFile)
temperatureProfile = [float(x) for x in config["Default"]["temperatureProfile"].split(",")]
#total heat given up by the liquid per second
flowRateKgPerSec = (float(config["Default"]["flow_rate"])/(60*1000))*float(config["Default"]["liquid_density"])
totalHeat = (float(config["Default"]["temp_in"])-float(config["Default"]["temp_out"]))*float(config["Default"]["c_liquid"])*flowRateKgPerSec
earTempTop = staveTop.getTemperatures("ear")[0]
if not args.one_face:
  earTempBottom = staveBottom.getTemperatures("ear")[0]

fractionHeat_segment0 = (temperatureProfile[1]-temperatureProfile[0])/(temperatureProfile[-1] - temperatureProfile[0])
fractionHeat_segment1 = (temperatureProfile[2]-temperatureProfile[1])/(temperatureProfile[-1] - temperatureProfile[0])
fractionHeat_segment2 = (temperatureProfile[3]-temperatureProfile[2])/(temperatureProfile[-1] - temperatureProfile[0])

logging.debug("fractionHeat_segment0 = {}".format(fractionHeat_segment0))
logging.debug("fractionHeat_segment1 = {}".format(fractionHeat_segment1))
logging.debug("fractionHeat_segment2 = {}".format(fractionHeat_segment2))

earHeat = (fractionHeat_segment0+fractionHeat_segment1 - 2*fractionHeat_segment2)*totalHeat/2
heatNextEar = (1.0 + 54.0/98)*fractionHeat_segment2*totalHeat/2
#liquid temperature between segments 0 and 1
liqTempAfterSeg0 = float(config["Default"]["temp_in"]) - temperatureProfile[1]*(float(config["Default"]["temp_in"])-float(config["Default"]["temp_out"]))

logging.debug("totalHeat = {}".format(totalHeat))
logging.debug("earTempTop = {}".format(earTempTop))
if not args.one_face:
  logging.debug("earTempBottom = {}".format(earTempBottom))
logging.debug("earHeat = {}".format(earHeat))
logging.debug("heatNextEar = {}".format(heatNextEar))
logging.debug("liqTempAfterSeg0 = {}".format(liqTempAfterSeg0))

#dT/dQ_region_segment as described in Documents/2020-09-09-EOS-Impedances.pdf
#importing the values from the config file
logging.debug("Loading the correction factors from the config file:")
dTdQ_large_0 = float(config["Default"]["dTdQ_large_0"]) #1.193
dTdQ_large_1 = float(config["Default"]["dTdQ_large_1"]) #0.716
dTdQ_small_0 = float(config["Default"]["dTdQ_small_0"]) #0.591
dTdQ_small_1 = float(config["Default"]["dTdQ_small_1"]) #0.251
dTdQ_nextEar = float(config["Default"]["dTdQ_nextEar"]) #1.152

logging.debug("dTdQ_large_0 = {}".format(dTdQ_large_0))
logging.debug("dTdQ_large_1 = {}".format(dTdQ_large_1))
logging.debug("dTdQ_small_0 = {}".format(dTdQ_small_0))
logging.debug("dTdQ_small_1 = {}".format(dTdQ_small_1))
logging.debug("dTdQ_nextEar = {}".format(dTdQ_nextEar))

#correcting the surface temperatures of the segments around the EoS region
staveTop.setTemperatureCorrection("large",0, earHeat*dTdQ_large_0)
staveTop.setTemperatureCorrection("large",1, earHeat*dTdQ_large_1)
staveTop.setTemperatureCorrection("small",0, earHeat*dTdQ_small_0)
staveTop.setTemperatureCorrection("small",1, earHeat*dTdQ_small_1)
if not args.one_face:
  staveBottom.setTemperatureCorrection("large",0, earHeat*dTdQ_large_0)
  staveBottom.setTemperatureCorrection("large",1, earHeat*dTdQ_large_1)
  staveBottom.setTemperatureCorrection("small",0, earHeat*dTdQ_small_0)
  staveBottom.setTemperatureCorrection("small",1, earHeat*dTdQ_small_1)

logging.debug("Temperature corrections for staveTop small regions: {}".format(str(staveTop.getTemperatureCorrections("small"))))
logging.debug("Temperature corrections for staveTop large regions: {}".format(str(staveTop.getTemperatureCorrections("large"))))
if not args.one_face:
  logging.debug("Temperature corrections for staveBottom small regions: {}".format(str(staveBottom.getTemperatureCorrections("small"))))
  logging.debug("Temperature corrections for staveBottom large regions: {}".format(str(staveBottom.getTemperatureCorrections("large"))))

#computing the impedance for the ear
earImpedanceTop = (liqTempAfterSeg0 - earTempTop - heatNextEar*dTdQ_nextEar)/earHeat
print("Z_earTop = {}".format(earImpedanceTop))

if not args.one_face:
  earImpedanceBottom = (liqTempAfterSeg0 - earTempBottom - heatNextEar*dTdQ_nextEar)/earHeat
  print("Z_earBottom = {}".format(earImpedanceBottom))

#WIP: print the impedance on the plot as well

#extracting the impedances
largeTop = staveTop.getImpedances("large", heatCorrection=True)
smallTop = staveTop.getImpedances("small")

if not(args.one_face):
  largeBottom = staveBottom.getImpedances("large", heatCorrection=True)
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
  f.write("\n")
  f.write("Z_earTop, {} \n".format(earImpedanceTop))
  if not args.one_face:
    f.write("Z_earBottom, {}".format(earImpedanceBottom))
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
    
  plt.plot([-1],[earImpedanceTop],marker='o', linestyle='', label="Z_earTop")
  if not args.one_face:
    plt.plot([-1],[earImpedanceBottom],marker='o', linestyle='', label="Z_earBottom")
  
  plt.xlabel("Region number")
  plt.ylabel("Thermal Impedance [K/W]")
  plt.title(outputFilename.split("/")[-1])
  if args.one_face:
    yrange = int(1+1.1*np.max([np.max(largeTop),np.max(smallTop),earImpedanceTop]))
  else:
    yrange = int(1+1.1*np.max([np.max(largeTop),np.max(largeBottom),np.max(smallTop),np.max(smallBottom),earImpedanceTop,earImpedanceBottom]))
  plt.xticks(np.arange(0, 28, 1.0))
  plt.yticks(np.arange(0, yrange, 0.5))
  plt.axis([-2.0,27.5,0,yrange])
  plt.grid()
  plt.legend(ncol=3)
  #ear impedances printed on the plot
  """
  ZearStr = "Z_earTop = {}".format(earImpedanceTop)
  if not args.one_face:
    ZearStr = ZearStr + "       Z_earBottom = {}".format(earImpedanceBottom)
  plt.text(0, 0.9*yrange, ZearStr, fontsize=10, bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
  """
  #code version printed on the plot
  plt.text(0, -0.13*yrange, "Code version: {} {}".format(gitHash, gitDate[:-6]), fontsize=10)
  plt.savefig(outputFilename + ".png")
  print("Outputing graphical output into a file: " + outputFilename + ".png")


if args.debug:
  plt.clf()
  plt.imshow(img_edges)
  plt.savefig("debug_output/edges.png")
