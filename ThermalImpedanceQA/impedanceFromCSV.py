#!/usr/bin/env python

'''
impedanceFromCSV.py

Authors: Lubos Vozdecky (Queen Mary, University of London)
About: This program takes a thermal image of an ATLAS Itk Stave Support in CSV format,
  finds the stave and computes the thermal impedances using data of the cooling fluid saved in parameters.cfg.
  The program outputs the result impedances into the /data folder as a CSV file.

Requires: pyROOT, Python 2.7, OpenCV

'''

import numpy as np
import cv2
import ROOT
import csv
import sys
import configparser
from matplotlib import pyplot as plt

def FindPoints(image, aproxYpos):
  """
  This function takes an input root stave image and finds all of the appropriate
  locations on the stave and creates a config file.
  """
  
  xPixels = image.shape[0]
  yPixels = image.shape[1]
  fltxPercentCutL=0.05
  fltxPercentCutR=0.023
  fltyPercentCut=0.20
  bol14ModCore = True
  
  # Make The Canny Image
  v = np.median(image)

  sigma = 0.33
  lower = int(max(0,(1-sigma)*v))
  upper = int(min(255,(1+sigma)*v))

  laplacian = cv2.Canny(np.uint8(image),lower,upper)
  image2 = laplacian
  """
  #Makes a Canny Image that can be checked
  histcanny = ROOT.TH2F("cannyplot","Canny Plot;xPixel;yPixel",xPixels,0,xPixels,yPixels,0,yPixels)
  orighist = ROOT.TH2F("originalplot","OriginalPlot;xPixel;yPixel",xPixels,0,xPixels,yPixels,0,yPixels)
  for i in range(xPixels):
    for j in range(yPixels):
      orighist.Fill(i,j,image[i][j])
      if image2[i][j] > 100:
        histcanny.Fill(i,j,image2[i][j])
      else:
        histcanny.Fill(i,j,0)
        image2[i][j] = 0 #Replace any small numbers with 0

  c2 = ROOT.TCanvas("c2")
  c2.cd()
  histcanny.Draw("colz")
  c2.Update()
 """

  # Find the Four Corners of the Pipe Area

 # c2.lines = []  #The lines stored in the canvas to show the buildup
  lineData = []  #This will be the four points
  HorData  = []  #The average y value of each horizontal line
  VertData = []  #The average x value of each vertical line
  ShortHorData = []

  #orighist.Draw("colz")
 
#------------------------------------------------------------------------------
  try:
    #Find Long Horiz Lines
    findLongLines = cv2.HoughLinesP(image2,rho = 1,theta = 1*np.pi/1000,threshold = 100,minLineLength = 200,maxLineGap = 175)

    findLongLines = findLongLines.flatten()
    nlines = len(findLongLines)/4
    findLongLines = np.split(findLongLines,nlines)

   
    LengthHoriz = nlines
    for line in range(int(LengthHoriz)):
      x1 = findLongLines[line][1]
      y1 = findLongLines[line][0]
      x2 = findLongLines[line][3]
      y2 = findLongLines[line][2]

                
      slope = abs((y2-y1)/(x2-x1+0.0001))
      intercept =y1-x1*(slope)

      HorData += [(y1+y2)/2]
      #lineObj = ROOT.TLine(x1,y1,x2,y2)
      #lineObj.SetLineColor(3)
      #lineObj.SetLineWidth(3)
      #c2.lines += [lineObj]
      #lineObj.Draw()

   
    HorData = np.sort(HorData)
    CentSep = np.amax(abs(HorData-aproxYpos)) 
    while CentSep > 50:
      lenHor = np.size(HorData)
      if abs(HorData[0]-aproxYpos)> 50:
        HorData = np.delete(HorData,0)
      else:
        HorData = np.delete(HorData,lenHor-1)
      CentSep = np.amax(abs(HorData-aproxYpos))

  
    lineData += [np.amax(HorData)]
    lineData += [np.amin(HorData)]

    
  except:
    print("Failed to Find Long Horizontal Lines. Using standard value for centered stave core")
    lineData += [265]
    lineData += [215]

#------------------------------------------------------------------------------   
#Find Short Vert Lines
  if bol14ModCore == False:
    StaveLength = 548.2 #13 module stave core length in pixels
    cutPercent= 0.05
  else:
    #StaveLength = 580 #~14 module stave core length (This is approximated from a Yale thermal image from Rec-000110)
    StaveLength = 630 #~14 module stave for the QMUL IR set-up
    #cutPercent= 0.02
    cutPercent = 0.20 # QMUL
  try:
    findShortLines = cv2.HoughLinesP(image2,rho = 1,theta = 1*np.pi/10000,threshold = 20,minLineLength = 10, maxLineGap = 5)

    findShortLines = findShortLines.flatten()
    nlines = len(findShortLines)/4
    findShortLines = np.split(findShortLines,nlines)

    LengthVert = nlines
    for line in range(int(LengthVert)):
      x1 = findShortLines[line][1]
      y1 = findShortLines[line][0]
      x2 = findShortLines[line][3]
      y2 = findShortLines[line][2]
      slope = abs((y2-y1)/(x2-x1+0.0001))
      intercept = y1-x1*(slope)

      if slope > 1: #Since this fit will find many short line segments we want to remove all horiztal lines
        VertData += [(x1+x2)/2]
        #lineObj = ROOT.TLine(x1,y1,x2,y2)
        #lineObj.SetLineColor(2)
        #lineObj.SetLineWidth(4)
        #c2.lines += [lineObj]
        #lineObj.Draw()

    VertData = np.sort(VertData)
    VertData = np.unique(VertData)
    
    #Check to confine to a stave core #Updated April 18th 2019 to make better x position choices
    VertSep = np.amax(VertData)-np.amin(VertData)
    skipConfining = False
    if VertSep < StaveLength- StaveLength*cutPercent:
      print("\n Unable to find all vertical lines around stave!!!! \n Trying to find best selection from found lines for the x position.\n")
      minPoint = np.amin(VertData)
      maxPoint = np.amax(VertData)

      #Found a line at the EOS stave end
      if minPoint + StaveLength < 640:
        #Find the best 
        while len(VertData) > 0:
          #Remove outside possible
          if maxPoint + StaveLength > 640:
            VertData = np.delete(VertData,-1)
            maxPoint = np.amax(VertData)
          #Get best point
          else:
            bestPoint = 0
            bestLineTemp = 0.0
            for Point in VertData:
              avgLineTemp=0.0
              for i in range(int(Point),int(Point+StaveLength)):
                avgLineTemp += abs(image[i][yPixels/2] - 20.)
              if avgLineTemp > bestLineTemp:
                bestPoint = Point
                bestLineTemp = avgLineTemp
            newX = bestPoint + int(StaveLength)
            VertData=[bestPoint,newX] 
            skipConfining = True
            break

      #Found a line at the stave end   
      elif maxPoint - StaveLength > 0:   
        #Find the best 
        while len(VertData) > 0:
          #Remove outside possible
          if minPoint - StaveLength > 640-StaveLength:
            VertData = np.delete(VertData,0)
            minPoint = np.amin(VertData)
          #Get best point
          else:
            bestPoint = 0
            bestLineTemp = 0.0
            for Point in VertData:
              avgLineTemp=0.0
              for i in range(int(Point-StaveLength),int(Point)):
                avgLineTemp += abs(image[i][yPixels/2] - 20.)
              if avgLineTemp > bestLineTemp:
                bestPoint = Point
                bestLineTemp = avgLineTemp
            newX = int(bestPoint-StaveLength)
            VertData=[newX,bestPoint]        
            skipConfining = True
            break

      #Found something else
      else:
        print("\n Found No Usable lines.")

    #Loop confining lines found on x to the stave. If this fails 
    while VertSep > StaveLength + StaveLength*cutPercent: #Finds the size of the stave core less than 5% its maximum 
      if skipConfining == True: break
      if len( VertData) == 2:
        VertSep = np.amax(VertData)-np.amin(VertData)
        if VertSep < StaveLength + StaveLength*cutPercent and VertSep > StaveLength - StaveLength*cutPercent: break
        else:
          print("Found poor separation of "+ str(VertSep))
          raise("Crud")
      print VertSep
      VertDataFrontRemoved = np.delete(VertData,0) 
      VertDataBackRemoved = np.delete(VertData,-1)
      VertSepFR = np.amax(VertDataFrontRemoved)-np.amin(VertDataFrontRemoved)     
      VertSepBR = np.amax(VertDataBackRemoved)-np.amin(VertDataBackRemoved)
      if VertSepFR > VertSepBR:
        VertSep = VertSepFR
        VertData = VertDataFrontRemoved
      else:
        VertSep = VertSepBR
        VertData = VertDataBackRemoved

    lineData += [np.amax(VertData)]
    lineData += [np.amin(VertData)]
   
  except:
    print("Failed to find Short Vertical Lines. Using different method assuming 13 module stave core length")
    VertData = [i for i in range(int(xPixels-StaveLength))]
    bestPoint = 0
    bestLineTemp = 0.0
    for Point in VertData:
      avgLineTemp=0.0
      for i in range(int(Point),int(Point+StaveLength)):
        avgLineTemp += abs(image[i][yPixels/2] - 20.)
      if avgLineTemp > bestLineTemp:
        bestPoint = Point
        bestLineTemp = avgLineTemp
    newX = bestPoint + int(StaveLength)        
    lineData+= [newX]
    lineData+= [bestPoint]  
#------------------------------------------------------------------------------
  """
  #Put the found area on the plot 
  for i in range(2):
    lineObj = ROOT.TLine(lineData[2],lineData[i],lineData[3],lineData[i])
    lineObj2 = ROOT.TLine(lineData[i+2],lineData[0],lineData[i+2],lineData[1])
    lineObj.SetLineColor(1)
    lineObj2.SetLineColor(1)
    lineObj.SetLineWidth(1)
    lineObj2.SetLineWidth(1)
    c2.lines += [lineObj]
    lineObj.Draw()
    c2.lines += [lineObj2]
    lineObj2.Draw()
   """



  x0=lineData[3]
  x1=lineData[2]
  y0=lineData[1]
  y1=lineData[0]
  Dx = x1-x0 
  Dy = y1-y0

  return [x0,y0,x1,y1]


if len(sys.argv) <= 1:
  print("Missing argument. First argument should be a .csv file")
  quit()

inputFile =  sys.argv[1]

#check if the suffix is .csv

if inputFile[-3:] != "csv":
  print("The input file should be a .csv file")
  quit()


#upload the CSV file 
image = []
with open(inputFile) as csvfile:
  reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    image.append(row)

img = np.array(image)

height = img.shape[0]
width = img.shape[1]

# split the image into two images
upper_img = img[0:int(height/2), 0:width].transpose()
lower_img = img[int(height/2):height, 0:width].transpose()


print("------------------------------------------------")

#look for the edges in both (uppe and lower) part of the images
#FindPoits function has to be given the approximate Y position of the stave core
points_upper = FindPoints(upper_img,200)
print("upper points: " + str(points_upper))


points_lower = FindPoints(lower_img,40)
print("lower points: " + str(points_lower))


print("------------------------------------------------")

#drawing the rectangles on the IR image
#deep copy of the image to draw on the image without affecting the original one that is used to compute the temperatures
img3 = np.copy(img)


#cv2.rectangle(img3,(points_upper[0],points_upper[1]),(points_upper[2],points_upper[3]),(0,0,0),thickness=1)
#cv2.rectangle(img3,(points_lower[0],points_lower[1]+int(height/2)),(points_lower[2],points_lower[3]+int(height/2)),(0,0,0),thickness=1)


#drawing the modules regions and calculating the average temperature
number_of_modules = 14
upper_stave_length = points_upper[2] - points_upper[0]
upper_stave_width = points_upper[3] - points_upper[1]
upper_module_length = upper_stave_length/(number_of_modules*1.0)  #need to convert to double to get rid of rounding error

topStave_largeUpperRegion_temp = []
topStave_largeLowerRegion_temp = []
bottomStave_largeUpperRegion_temp = []
bottomStave_largeLowerRegion_temp = []

topStave_smallUpperRegion_temp = []
topStave_smallLowerRegion_temp = []
bottomStave_smallUpperRegion_temp = []
bottomStave_smallLowerRegion_temp = []


#calculating number of pixels (vertically) for the small region (8.0mm wide)
pixelsSmallRegion = int(round((4.0*upper_stave_width)/115))

if pixelsSmallRegion == 0:
  pixelsSmallRegion=1




for i in range(1,number_of_modules+1):
  #cv2.rectangle(img3,(points_upper[0]+int((i-1)*upper_module_length),points_upper[1]),(points_upper[0]+int(i*upper_module_length),points_upper[1]+upper_stave_width/2),(0,0,0),thickness=1)
  #cv2.rectangle(img3,(points_upper[0]+int((i-1)*upper_module_length),points_upper[1]+upper_stave_width/2),(points_upper[0]+int(i*upper_module_length),points_upper[3]),(0,0,0),thickness=1)

  cv2.rectangle(img3,(points_upper[0]+int((i-1)*upper_module_length),points_upper[1]+ int(0.2826*upper_stave_width)-pixelsSmallRegion),(points_upper[0]+int(i*upper_module_length),points_upper[1]+ int(0.2826*upper_stave_width)+pixelsSmallRegion),(0,0,0),thickness=1)
  cv2.rectangle(img3,(points_upper[0]+int((i-1)*upper_module_length),points_upper[1]+ int(upper_stave_width/2+0.2826*upper_stave_width)-pixelsSmallRegion),(points_upper[0]+int(i*upper_module_length),points_upper[1]+ int(upper_stave_width/2+0.2826*upper_stave_width)+pixelsSmallRegion),(0,0,0),thickness=1)
  
  top_y = points_upper[1]
  bottom_y = points_upper[1]+upper_stave_width/2
  left_x = points_upper[0]+int((i-1)*upper_module_length)
  right_x = points_upper[0]+int(i*upper_module_length)

  #print("[" + str(left_x) + ", " + str(top_y) + "]" +  " ..... [" +  str(right_x) + ", " + str(bottom_y) + "]")
  
  crop = img[top_y:bottom_y,left_x:right_x]
  topStave_largeUpperRegion_temp.append(np.mean(crop))
  
  
  top_y = points_upper[1]+upper_stave_width/2
  bottom_y = points_upper[3]
  left_x = points_upper[0]+int((i-1)*upper_module_length)
  right_x = points_upper[0]+int(i*upper_module_length)

  crop = img[top_y:bottom_y,left_x:right_x]
  topStave_largeLowerRegion_temp.append(np.mean(crop))
  
  
  
  top_y = points_upper[1]+ int(0.2826*upper_stave_width)-pixelsSmallRegion
  bottom_y = points_upper[1]+ int(0.2826*upper_stave_width)+pixelsSmallRegion
  left_x = points_upper[0]+int((i-1)*upper_module_length)
  right_x = points_upper[0]+int(i*upper_module_length)
  
  crop = img[top_y:bottom_y,left_x:right_x]
  topStave_smallUpperRegion_temp.append(np.mean(crop))
  
  
  top_y = points_upper[1]+ int(upper_stave_width/2+0.2826*upper_stave_width)-pixelsSmallRegion
  bottom_y = points_upper[1]+ int(upper_stave_width/2+0.2826*upper_stave_width)+pixelsSmallRegion
  left_x = points_upper[0]+int((i-1)*upper_module_length)
  right_x = points_upper[0]+int(i*upper_module_length)
  
  crop = img[top_y:bottom_y,left_x:right_x]
  topStave_smallLowerRegion_temp.append(np.mean(crop))


lower_stave_length = points_lower[2] - points_lower[0]
lower_stave_width = points_lower[3] - points_lower[1]
lower_module_length = lower_stave_length/(number_of_modules*1.0) #need to convert to double to get rid of rounding error


for i in range(1,number_of_modules+1):
  #cv2.rectangle(img3,(points_lower[0]+int((i-1)*lower_module_length),points_lower[1]+int(height/2)),(points_lower[0]+int(i*lower_module_length),points_lower[1]+int(height/2)+lower_stave_width/2),(0,0,0),thickness=1)
  #cv2.rectangle(img3,(points_lower[0]+int((i-1)*lower_module_length),points_lower[1]+int(height/2)+lower_stave_width/2),(points_lower[0]+int(i*lower_module_length),points_lower[3]+int(height/2)),(0,0,0),thickness=1)

  cv2.rectangle(img3,(points_lower[0]+int((i-1)*lower_module_length),points_lower[1]+int(height/2)+int(0.2826*lower_stave_width)-pixelsSmallRegion),(points_lower[0]+int(i*lower_module_length),points_lower[1]+int(height/2)+int(0.2826*lower_stave_width)+pixelsSmallRegion),(0,0,0),thickness=1)
  cv2.rectangle(img3,(points_lower[0]+int((i-1)*lower_module_length),points_lower[1]+int(height/2)+lower_stave_width/2+int(0.2826*lower_stave_width)-pixelsSmallRegion),(points_lower[0]+int(i*lower_module_length),points_lower[1]+int(height/2)+lower_stave_width/2+int(0.2826*lower_stave_width)+pixelsSmallRegion),(0,0,0),thickness=1)

  top_y = points_lower[1]+int(height/2)
  bottom_y = points_lower[1]+lower_stave_width/2+int(height/2)
  left_x = points_lower[0]+int((i-1)*lower_module_length)
  right_x = points_lower[0]+int(i*lower_module_length)
  
  crop = img[top_y:bottom_y,left_x:right_x]
  bottomStave_largeUpperRegion_temp.append(np.mean(crop))
  

  top_y = points_lower[1]+lower_stave_width/2+int(height/2)
  bottom_y = points_lower[3]+int(height/2)
  left_x = points_lower[0]+int((i-1)*lower_module_length)
  right_x = points_lower[0]+int(i*lower_module_length)

  #print("[" + str(left_x) + ", " + str(top_y) + "]" +  " ..... [" +  str(right_x) + ", " + str(bottom_y) + "]")
  crop = img[top_y:bottom_y,left_x:right_x]
  bottomStave_largeLowerRegion_temp.append(np.mean(crop))
  
  top_y = points_lower[1]+int(height/2)+int(0.2826*lower_stave_width)-pixelsSmallRegion
  bottom_y = points_lower[1]+int(height/2)+int(0.2826*lower_stave_width)+pixelsSmallRegion
  left_x = points_lower[0]+int((i-1)*lower_module_length)
  right_x = points_lower[0]+int(i*lower_module_length)
  
  crop = img[top_y:bottom_y,left_x:right_x]
  bottomStave_smallUpperRegion_temp.append(np.mean(crop))
  
  
  top_y = points_lower[1]+int(height/2)+lower_stave_width/2+int(0.2826*lower_stave_width)-pixelsSmallRegion
  bottom_y = points_lower[1]+int(height/2)+lower_stave_width/2+int(0.2826*lower_stave_width)+pixelsSmallRegion
  left_x = points_lower[0]+int((i-1)*lower_module_length)
  right_x = points_lower[0]+int(i*lower_module_length)
  
  crop = img[top_y:bottom_y,left_x:right_x]
  bottomStave_smallLowerRegion_temp.append(np.mean(crop))
 

print("------------------------------------------------")

#importing data from parameters.cfg
filename = "parameters.cfg"
print("Importing variables from config file " + filename + ":")
config = configparser.ConfigParser()
config.read(filename)

Tin = float(config["Default"]["temp_in"])
Tout = float(config["Default"]["temp_out"])
Cliq = float(config["Default"]["c_liquid"])
Zcut = float(config["Default"]["z_cut"])
FR = float(config["Default"]["flow_rate"])

#print the imported variables
for variable in config.items("Default"):
  print(variable[0] + " = " + variable[1]) 
print("------------------------------------------------")


#temperature profile of the liquid comes from the finite element analysis of the stave
temperatureProfile = [0,0.062,0.114,0.152,0.188,0.224,0.26,0.295,0.33,0.364,0.398,0.432,0.466,0.499,0.533,0.568,0.603,0.637,0.672,0.706,0.74,0.774,0.807,0.841,0.873,0.906,0.937,0.969,1]


#scale scale it up
Tliquid = map(lambda x:x*((Tout-Tin)/temperatureProfile[-1]), temperatureProfile)

#shift it to match Tin
Tliquid = map(lambda x:x+(Tin-Tliquid[0]), Tliquid)

"""
#linear extrapolation of the temperature of the liquid
Tliquid = Tin - ((Tin-Tout)/28.0)*np.linspace(0,28,29)
"""

heatCapacity = Cliq

thermalImpedance_topLargeRegion = []
thermalImpedance_bottomLargeRegion = []
thermalImpedance_topSmallRegion = []
thermalImpedance_bottomSmallRegion = []


flowRate_kgPerSec = FR/60


#impedance for the large region
for i in range(0,14):
  # facor 0.5 for considering only one side
  heat = (Tliquid[i] - Tliquid[i+1])*heatCapacity*0.5*flowRate_kgPerSec
  averageTempDiff_upperFace  = (Tliquid[i] + Tliquid[i+1])/2 - topStave_largeUpperRegion_temp[i]
  #flow of cooling liquid is flipped for the bottom face as the thermal image is flipped => use (27-i) instead of i
  averageTempDiff_lowerFace  = (Tliquid[(27-i)] + Tliquid[(27-i)+1])/2 - bottomStave_largeUpperRegion_temp[i]
  thermalImpedance_topLargeRegion.append(averageTempDiff_upperFace/heat)
  thermalImpedance_bottomLargeRegion.append(averageTempDiff_lowerFace/heat)

for i in range(14,28):
  heat = (Tliquid[i] - Tliquid[i+1])*heatCapacity*0.5*flowRate_kgPerSec
  averageTempDiff_upperFace  = (Tliquid[i] + Tliquid[i+1])/2 - topStave_largeLowerRegion_temp[i-14]
  #flow of cooling liquid is flipped for the bottom face as the thermal image is flipped => use (27-i) instead of i
  averageTempDiff_lowerFace = (Tliquid[(27-i)] + Tliquid[(27-i)+1])/2 - bottomStave_largeLowerRegion_temp[i-14]
  thermalImpedance_topLargeRegion.append(averageTempDiff_upperFace/heat)
  thermalImpedance_bottomLargeRegion.append(averageTempDiff_lowerFace/heat)
  

#impedance for the small region
for i in range(0,14):
  # facor 0.5 for considering only one side
  heat = (Tliquid[i] - Tliquid[i+1])*heatCapacity*0.5*flowRate_kgPerSec
  averageTempDiff_upperFace  = (Tliquid[i] + Tliquid[i+1])/2 - topStave_smallUpperRegion_temp[i]
  #flow of cooling liquid is flipped for the bottom face as the thermal image is flipped => use (27-i) instead of i
  averageTempDiff_lowerFace  = (Tliquid[(27-i)] + Tliquid[(27-i)+1])/2 - bottomStave_smallUpperRegion_temp[i]
  thermalImpedance_topSmallRegion.append(averageTempDiff_upperFace/heat)
  thermalImpedance_bottomSmallRegion.append(averageTempDiff_lowerFace/heat)

for i in range(14,28):
  heat = (Tliquid[i] - Tliquid[i+1])*heatCapacity*0.5*flowRate_kgPerSec
  averageTempDiff_upperFace  = (Tliquid[i] + Tliquid[i+1])/2 - topStave_smallLowerRegion_temp[i-14]
  #flow of cooling liquid is flipped for the bottom face as the thermal image is flipped => use (27-i) instead of i
  averageTempDiff_lowerFace = (Tliquid[(27-i)] + Tliquid[(27-i)+1])/2 - bottomStave_smallLowerRegion_temp[i-14]
  thermalImpedance_topSmallRegion.append(averageTempDiff_upperFace/heat)
  thermalImpedance_bottomSmallRegion.append(averageTempDiff_lowerFace/heat)




"""
print("Upper face:")
print(thermalImpedance_topLargeRegion)
print("Lower face:")
print(thermalImpedance_bottomLargeRegion)
"""

"""
print("Large region, Upper Face:")
print("# \t Impedance \t QA")
for i in range(0,28):
  print(str(i) + "\t" + "%0.3f" % thermalImpedance_topLargeRegion[i] + "\t\t" + ("OK" if thermalImpedance_topLargeRegion[i] <= Zcut else "FAIL"))

print("-------------------------")

print("Large region, Lower Face:")
print("# \t Impedance \t QA")
for i in range(0,28):
  print(str(i) + "\t" + "%0.3f" % thermalImpedance_bottomLargeRegion[i] + "\t\t" + ("OK" if thermalImpedance_bottomLargeRegion[i] <= Zcut else "FAIL"))


print("-------------------------")
print("-------------------------")
print("-------------------------")

print("Small region, Upper Face:")
print("# \t Impedance \t QA")
for i in range(0,28):
  print(str(i) + "\t" + "%0.3f" % thermalImpedance_topSmallRegion[i] + "\t\t" + ("OK" if thermalImpedance_topSmallRegion[i] <= Zcut else "FAIL"))

print("-------------------------")


print("Small region, Lower Face:")
print("# \t Impedance \t QA")
for i in range(0,28):
  print(str(i) + "\t" + "%0.3f" % thermalImpedance_bottomSmallRegion[i] + "\t\t" + ("OK" if thermalImpedance_bottomSmallRegion[i] <= Zcut else "FAIL"))

print("-------------------------")
print("-------------------------")
print("-------------------------")
print("topStave_smallUpperRegion_temp = " + str(topStave_smallUpperRegion_temp))
print("-------------------------")
print("topStave_smallLowerRegion_temp = " + str(topStave_smallLowerRegion_temp))
print("-------------------------")
print("bottomStave_smallUpperRegion_temp = " + str(bottomStave_smallUpperRegion_temp))
print("-------------------------")
print("bottomStave_smallLowerRegion_temp = " + str(bottomStave_smallLowerRegion_temp))
print("-------------------------")
print("Tliquid = " + str(Tliquid))
print("-------------------------")
"""


#plt.imshow(upper_img)
#plt.show()
#plt.imshow(lower_img)
#plt.show()

plt.imshow(img3)
plt.show()


#plt.imshow(crop_img, cmap='hot', interpolation='nearest')
#plt.show()

outputFilename = "output/" + inputFile.split("/")[-1][:-4] + "_IMPEDANCES.csv"
print("Outputing data into a file: " + outputFilename)
with open(outputFilename, "w+") as f:
  f.write("#, topLargeRegion, bottomLargeRegion, topSmallRegion, bottomSmallRegion \n")
  for i in range(0,28):
    f.write(str(i)+", "+str(thermalImpedance_topLargeRegion[i])+", "+str(thermalImpedance_bottomLargeRegion[i])+", "+str(thermalImpedance_topSmallRegion[i])+", "+str(thermalImpedance_bottomSmallRegion[i]) + "\n")

f.close()

""" 
  #Print all of the Fit Lines
  c2.Update()
  c2.Print(outdir+"/AllFoundLines.pdf")
  c2.Print(outdir+"/AllFoundLines.root")

  #Get Corrected Zoomed Figure
  DxcutL = int(fltxPercentCutL*Dx)
  DxcutR = int(fltxPercentCutR*Dx)
  Dycut = int(fltyPercentCut*Dy)

  x0Cut = x0+DxcutL
  x1Cut = x1-DxcutR
  y0Cut = y0+Dycut
  y1Cut = y1-Dycut
 
  c2.Close()

  #Make the output file
  Output = open(strOutputFile,"w")
  Output.write("#\n# frame parameters used in frameanal.py\n#\n")

  #Find the EndofStaveCard
  avgStaveTemp = 0.0
  avgAboveTemp = 0.0
  avgBelowTemp = 0.0
  stavePixels = (x1Cut-x0Cut)*(y1Cut-y0Cut)
  EOSCardPixels = 300 
  for i in range(xPixels):
    for j in range(yPixels):
      if i > x0Cut and i < x1Cut:
        if j > y0Cut and j < y1Cut:
          avgStaveTemp = avgStaveTemp + image[i][j]
      if i > x0Cut and i < (x0Cut+30):
        if j > y1 and j < (y1+10):
          avgAboveTemp += image[i][j]
        if j > y0 - 10 and j < y0:
          avgBelowTemp += image[i][j] 
      
  avgStaveTemp = avgStaveTemp/stavePixels
  avgAboveTemp = avgAboveTemp/EOSCardPixels
  avgBelowTemp = avgBelowTemp/EOSCardPixels

  #Get Which Side from the Plot
  if avgStaveTemp > 20: 
    if avgAboveTemp > avgBelowTemp:
      intStaveSideL = 0
    else:
      intStaveSideL = 1
  else:
    if avgAboveTemp < avgBelowTemp:
      intStaveSideL = 0
    else:
      intStaveSideL = 1

  if intStaveSideL == 0:
    y0EOS = y0
    y1EOS = y1 + int((y1-y0)*0.4)
  else: 
    y0EOS = y0 - int((y1-y0)*0.4)
    y1EOS = y1 

  if intStaveSideL ==1:
    yorig = [y0EOS,y0Cut,y1EOS,y1Cut]
    y0EOS = abs(yorig[2] - yPixels) 
    y0Cut = abs(yorig[3] - yPixels)
    y1EOS = abs(yorig[0] - yPixels)
    y1Cut = abs(yorig[1] - yPixels)

  #Put in the stave parameters
  Output.write("StavePixelX0 "+str(x0)+"\n")
  Output.write("StavePixelY0 "+str(y0EOS)+"\n")
  Output.write("StavePixelX1 "+str(x1)+"\n")
  Output.write("StavePixelY1 "+str(y1EOS)+"\n")

  #Put in the pipe area parameters
  Output.write("PipePixelX0 "+str(x0Cut)+"\n")
  Output.write("PipePixelY0 "+str(y0Cut)+"\n")
  Output.write("PipePixelX1 "+str(x1Cut)+"\n")
  Output.write("PipePixelY1 "+str(y1Cut)+"\n")

  #Put in the Other Constants
  Output.write("CMperPixel 0.23272727\n")
  Output.write("StaveSideL "+str(intStaveSideL)+"\n")    
  
  #State Whether it is hot or cold
  if avgStaveTemp > 25.:
    intLowTValue = 0
    #Assuming Run at 50C
    intMaxPlotTemp  = 50
    intMinPlotTemp  = 10
    intMaxStaveTemp = 50
    intMinStaveTemp = 25
    intMaxPipeTemp  = 50
    intMinPipeTemp  = 25

  elif avgStaveTemp <15.:
    intLowTValue = 1
    #Assuming Run at -55C
    intMaxPlotTemp  = 20
    intMinPlotTemp  = -40
    intMaxStaveTemp = -10
    intMinStaveTemp = -40
    intMaxPipeTemp  = -20
    intMinPipeTemp  = -40
  else:
    intLowTValue = 2 #This will cause the fits to not happen
    #Assuming Just Running Around Room Temp
    intMaxPlotTemp  = 25
    intMinPlotTemp  = 15
    intMaxStaveTemp = 25
    intMinStaveTemp = 15
    intMaxPipeTemp  = 15
    intMinPipeTemp  = 25

  Output.write("LiquidTLow "+str(intLowTValue)+"\n")

  #From Temp Values
  Output.write("FrameTmax "+str(intMaxPlotTemp)+"\n")
  Output.write("FrameTmin "+str(intMinPlotTemp)+"\n")
  Output.write("StaveTmax "+str(intMaxStaveTemp)+"\n")
  Output.write("StaveTmin "+str(intMinStaveTemp)+"\n")
  Output.write("PipeTmax "+str(intMaxPipeTemp)+"\n")
  Output.write("PipeTmin "+str(intMinPipeTemp)+"\n")

  Output.close()
  return [x0Cut,x1Cut,y0Cut,y1Cut,intStaveSideL,intLowTValue]

#Below is a short program that will get all of the fit data out of a set of frames

if __name__ == '__main__':
  x0hist = ROOT.TH1F("x0hist","X0 values found;xPixel;N Entries",10,89,99) 
  x1hist = ROOT.TH1F("x1hist","X1 values found;xPixel;N Entries",10,598,608) 
  y0hist = ROOT.TH1F("y0hist","Y0 values found;yPixel;N Entries",10,223,233) 
  y1hist = ROOT.TH1F("y1hist","Y1 values found;yPixel;N Entries",10,257,267) 
  staveSideHist = ROOT.TH1F("staveSideHist","Side Found;LSide?;N Entries",2,0,1) 
  TempHist = ROOT.TH1F("TempHist","AvgTemperature on Stave;Temperature [#circC];N Entries",200,-50,50) 
  for i in range(200):
    Output = FindPoints("roo/frame_"+str(i)+".root","configure")
    x0hist.Fill(Output[0])
    x1hist.Fill(Output[1])
    y0hist.Fill(Output[2])
    y1hist.Fill(Output[3])
    staveSideHist.Fill(Output[4])
    TempHist.Fill(Output[5])
    print(str(i)+" " +str(Output[0])+" "+str(Output[1])+" "+str(Output[2])+" "+str(Output[3])+" "+str(Output[4])+" "+str(Output[5]))

  c1= ROOT.TCanvas("c1")
  c1.cd()    
  x0hist.Draw()
  c1.Print("x0hist.png")
  x1hist.Draw()
  c1.Print("x1hist.png")
  y0hist.Draw()
  c1.Print("y0hist.png")
  y1hist.Draw()
  c1.Print("y1hist.png")
"""