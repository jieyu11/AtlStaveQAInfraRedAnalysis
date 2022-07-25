# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:16:24 2020

@author: Lubos (vozdeckyl@gmail.com)

Script that converts the FEA Abaqus output into a CSV file.
A file with settings named "settings.txt" is required.
"""

settings = []

with open("settings.txt", "r") as file:
    for line in file:
        if len(line) < 2:
            break
        settings.append(line.split("\t")[1])
        

rptFileStartLine = int(settings[0])
ambientTemperature = float(settings[1])
paddingX = int(settings[2])
paddingY = int(settings[3])
latticeXperiod = int(settings[4])
latticeYperiod = int(settings[5])


tempLinesSplit = []

with open("abaqus.rpt", "r") as file:
    lines = file.readlines()[19:]
    for line in lines:
        split = [s for s in line.split(" ") if s!=""]
        if len(split) <= 1:
            break
        tempLinesSplit.append([int(split[0]), float(split[1])])


export = False

coordLines = []

with open("Core20LS_2_var2.inp", "r") as file:
    for line in file:
        if "name=PolyShell3dPlanar" in line:
            export = True
            
        if export and "*Element" in line:
            break
        
        if export:
            coordLines.append(line)
            

coordLines = coordLines[2:]

coordLinesSplit = []

for line in coordLines:
    split = line.replace(" ","").split(",")
    coordLinesSplit.append([int(split[0]),int(float(split[1])/latticeXperiod),int(float(split[2])/latticeYperiod)])

import numpy as np

coordinates = np.array(coordLinesSplit)
temperatures = np.array(tempLinesSplit)


x_max = np.max(coordinates, axis=0)[1]
y_max = np.max(coordinates, axis=0)[2]

image = ambientTemperature*np.ones([x_max+1,y_max+1])

for coord in coordinates:
    nodeNum = coord[0]
    nodeX = coord[1]
    nodeY = coord[2]
    nodeTemp = temperatures[nodeNum-1][1]
    
    if not temperatures[nodeNum-1][0] == nodeNum:
        print("ERROR!!!!!")
        break
    
    image[nodeX,nodeY] = nodeTemp
    

image = np.rot90(image,1,axes=(0,1))

import cv2

imageResized = cv2.resize(image, (5*image.shape[1],2*image.shape[0]), interpolation = cv2.INTER_LINEAR)

imgPadded = np.pad(imageResized,((paddingY,paddingY),(paddingX,paddingX)), 'constant',constant_values=((ambientTemperature,ambientTemperature),(ambientTemperature,ambientTemperature)))

np.set_printoptions(suppress=True)
np.savetxt('output.csv', imgPadded, delimiter=',',fmt='%2.6f')

