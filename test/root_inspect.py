#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:56:46 2020

@author: vozdecky
"""
import os
import numpy as np
import ROOT
from matplotlib import pyplot as plt





xPixels = 640
yPixels = 480
strImageFile = "../roo/frame_average.root"
try:
    #Load in the file
    imageFile = ROOT.TFile(strImageFile,"read")
except:
    print("Failed to Load ImageFile")
    quit()
Tree = imageFile.Get("atree")
_temperature = np.zeros(1,dtype=float)
Tree.SetBranchAddress("temperature",_temperature)

#Load the image from TTree
image = np.full((xPixels,yPixels),-999.) #A tree full of -999 used as a placeholder

bolAvgFrame = False
if "frame_average.root" in strImageFile:
    bolAvgFrame = True
  
if bolAvgFrame == True:
    for i in range(xPixels):
        for j in range(yPixels):
            Tree.GetEntry(i*yPixels + j)      #Reading from an average frame
            image[i][j] = _temperature[0] 

else:
    for i in range(xPixels):
        for j in range(yPixels):
            Tree.GetEntry(j*xPixels + i) #Reading from a single frame
            image[i][j] = _temperature[0]
            

np.savetxt("A655-000277-312_16_15_05_892_RAW.csv", image, delimiter=",")