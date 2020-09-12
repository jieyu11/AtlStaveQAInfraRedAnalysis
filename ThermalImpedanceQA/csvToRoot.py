#!/usr/bin/env python

import numpy as np
import os
import csv
import logging
import argparse
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument("path", help="The path to the input CSV file")
parser.add_argument("-d","--debug", help="Runs the code in debug mode", action="store_true")
args = parser.parse_args()

inputFile = args.path

#check if the suffix is .csv
if inputFile[-3:] != "csv":
  print("The input file should be a .csv file")
  quit()

#set up the debugging log
if args.debug:
  logging.basicConfig(filename='debug_output/csvToRoot_debug.log',level=logging.DEBUG)

#fetch the CSV file
logging.debug("Opening the CSV file")
imgList = []
with open(inputFile) as csvfile:
  reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    imgList.append(row)
image = np.array(imgList)

file = ROOT.TFile( "test.root", "recreate")

temperature = np.empty((1), dtype="float64")
xpos = np.empty((1), dtype="int64")
ypos = np.empty((1), dtype="int64")


atree = ROOT.TTree("atree", "a tree of temperature data")
atree.Branch('temperature', temperature, 'temperature/D')
atree.Branch('xpos', xpos, 'xpos/I')
atree.Branch('ypos', ypos, 'ypos/I')

for i in range(0,image.shape[0]):
    for j in range(0,image.shape[1]):
        temperature[0] = image[i][j]
        xpos[0] = i
        ypos[0] = j
        atree.Fill()

file.Write()
