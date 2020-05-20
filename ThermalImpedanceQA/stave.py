#!/usr/bin/env python

import logging
import os
import cv2
import configparser
import numpy as np
import math
from matplotlib import pyplot as plt

class Stave:
  
  def __init__(self,globalImg, configFile, args):
    self.__globalImg = globalImg
    self.__staveFound = False
    self.__xLeft = 0
    self.__xRight = 0
    self.__yTop = 0
    self.__yBottom = 0
    self.__regions = {} # dictionary
    self.__args = args
    self.__staveLength = 630 #~14 module stave for the QMUL IR set-up
    self.__cutPercentLength = 0.05 # 5% tolerance on the stave length
    #parameters for the edge finding algorithms:
    self.__minHorizontalLineLength = 50
    self.__minVerticalLineLength = 10
    self.__maxHorizontalLineGap = 10
    self.__maxVerticalLineGap = 5
    

    #importing data from the config file
    if not os.path.isfile(configFile):
      raise Exception("Could not find config file '" + configFile + "'")
	
    logging.debug("Importing variables from config file " + configFile + ":")
    config = configparser.ConfigParser()
    config.read(configFile)

    self.__Tin = float(config["Default"]["temp_in"])
    self.__Tout = float(config["Default"]["temp_out"])
    self.__heatCapacity = float(config["Default"]["c_liquid"])
    self.__Zcut = float(config["Default"]["z_cut"])
    self.__FR = float(config["Default"]["flow_rate"])
    self.__region_relaxation = bool(int(config["Default"]["region_relaxation"]))
    
    #print the imported variables
    for variable in config.items("Default"):
      logging.debug(variable[0] + " = " + variable[1])
      
    #default temperature profile of the cooling liquid (linear extraploation)
    self.__temperatureProfile = [x/28.0 for x in list(range(0,29))]
    
  def FindStaveWithin(self, relXMin, relXMax, relYMin, relYMax):
    if relXMin > relXMax or relYMin > relYMax:
      raise Exception("FindStaveWithin(): minimal values are larger than the maximum.")
    
    if self.__staveFound:
      raise Exception("Stave was already found. Cannot call FindStaveWithin again.")
    
    imageHeight = self.__globalImg.shape[0]
    imageWidth = self.__globalImg.shape[1]
    
    xMin = int(relXMin * imageWidth)
    xMax = int(relXMax * imageWidth)
    yMin = int(relYMin * imageHeight)
    yMax = int(relYMax * imageHeight)
    
    imgOfInterest = self.__globalImg[yMin:yMax,xMin:xMax]
    
    xPixels = imgOfInterest.shape[0]
    yPixels = imgOfInterest.shape[1]
    
    # Make The Canny Image
    v = np.median(imgOfInterest)
    #sigma = 0.33
    sigma = 0.05
    
    print("max: " + str(np.max(imgOfInterest)))
    print("min: " + str(np.min(imgOfInterest)))
    
    lower = int(max(0,(1-sigma)*v))
    upper = int(min(255,(1+sigma)*v))
    
    cannyImage = cv2.Canny(np.uint8(imgOfInterest),10,10)
    
    plt.imshow(cannyImage)
    plt.show()

  
    #sva the canny image into the debug folder
    if self.__args.debug:
      plt.figure(figsize=(12,12))
      plt.imshow(cannyImage)
      if os.path.isfile("debug_output/canny.png"):
        plt.savefig("debug_output/canny2.png")
      else:
        plt.savefig("debug_output/canny.png")
        plt.clf()
    
    
    # Find the Four Corners of the Pipe Area:
    
    lineData = []  #This will be the four points
    HorData  = []  #The average y value of each horizontal line
    VertData = []  #The average x value of each vertical line
    ShortHorData = []
    
    try:
      #Find Long Horiz Lines
      findLongLines = cv2.HoughLinesP(cannyImage,rho = 1,theta = 1*np.pi/1000,threshold = 20,minLineLength = self.__minHorizontalLineLength,maxLineGap = self.__maxHorizontalLineGap)
      findLongLines = findLongLines.flatten()
      nlines = len(findLongLines)/4
      findLongLines = np.split(findLongLines,nlines)
      
      print("longlines count = " + str(len(findLongLines)))
      
      LengthHoriz = nlines
      for line in range(int(LengthHoriz)):
        y1 = findLongLines[line][1]
        x1 = findLongLines[line][0]
        y2 = findLongLines[line][3]
        x2 = findLongLines[line][2]
                  
        slope = abs((y2-y1)/(x2-x1+0.0001))
        
        if abs(slope) > 0.1:
          continue
  
        intercept =y1-x1*(slope)
        HorData += [(y1+y2)/2]

      aproxYpos = np.mean(HorData)
      
      print("aproxYpos = " + str(aproxYpos))
      print("HorData = " + str(HorData))
      
      upperEdge = np.max(np.extract(HorData<=aproxYpos,HorData)[0])      
      lowerEdge = np.min(np.extract(HorData>aproxYpos,HorData)[0])

      lineData += [lowerEdge]
      lineData += [upperEdge]
      
      print("lowerEdge = " + str(lowerEdge))
      print("upperEdge = " + str(upperEdge))
    
    except:
      logging.error("Failed to Find Long Horizontal Lines")
      return
    
    
    
    StaveLength = self.__staveLength 
    cutPercent = self.__cutPercentLength
    
    #Find Short Vert Lines
    try:
      findShortLines = cv2.HoughLinesP(cannyImage,rho = 1,theta = 1*np.pi/10000,threshold = 20,minLineLength = self.__minVerticalLineLength, maxLineGap = self.__maxVerticalLineGap)

      findShortLines = findShortLines.flatten()
      nlines = len(findShortLines)/4
      findShortLines = np.split(findShortLines,nlines)

      LengthVert = nlines
      for line in range(int(LengthVert)):
        y1 = findShortLines[line][1]
        x1 = findShortLines[line][0]
        y2 = findShortLines[line][3]
        x2 = findShortLines[line][2]
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
        logging.debug("\n Unable to find all vertical lines around stave!!!! \n Trying to find best selection from found lines for the x position.\n")
        minPoint = np.amin(VertData)
        maxPoint = np.amax(VertData)

        #Found a line at the EOS stave end
        if minPoint + StaveLength < self.__globalImg.shape[1]:
          #Find the best 
          while len(VertData) > 0:
            #Remove outside possible
            if maxPoint + StaveLength > self.__globalImg.shape[1]:
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
            if minPoint - StaveLength > self.__globalImg.shape[1]-StaveLength:
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
          logging.debug("\n Found No Usable lines.")

      #Loop confining lines found on x to the stave. If this fails 
      while VertSep > StaveLength + StaveLength*cutPercent: #Finds the size of the stave core less than 5% its maximum 
        if skipConfining == True: break
        if len( VertData) == 2:
          VertSep = np.amax(VertData)-np.amin(VertData)
          if VertSep < StaveLength + StaveLength*cutPercent and VertSep > StaveLength - StaveLength*cutPercent: break
          else:
            logging.debug("Found poor separation of "+ str(VertSep))
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
      logging.error("Failed to find Short Vertical Lines.")
      return
    
    self.__xLeft = lineData[3] + xMin
    self.__xRight = lineData[2] + xMin
    self.__yTop = lineData[1] + yMin
    self.__yBottom = lineData[0] +yMin
    
    self.__width = self.__yBottom - self.__yTop
    self.__length = self.__xRight - self.__xLeft

    self.__staveFound = True
    self.__staveImg = self.__globalImg[self.__xLeft:self.__xRight,self.__yTop:self.__yBottom]


  def ScaleImage(self, scale):
    width = int(self.__globalImg.shape[1] * scale)
    height = int(self.__globalImg.shape[0] * scale)
    dim = (width, height)
    resized = cv2.resize(self.__globalImg, dim, interpolation = cv2.INTER_LINEAR)
    
    #scale up the parameters
    self.__staveLength = self.__staveLength * scale
    #parameters for the edge finding algorithms:
    self.__minHorizontalLineLength = self.__minHorizontalLineLength * scale
    self.__minVerticalLineLength = self.__minVerticalLineLength * scale
    self.__maxHorizontalLineGap = self.__maxHorizontalLineGap * scale
    self.__maxVerticalLineGap = self.__maxVerticalLineGap * scale
    
    self.__globalImg = resized
  
  
  def setTemperatureProfile(self,newTemperatureProfile):
    self.__temperatureProfile = newTemperatureProfile

  def AddRegion(self,xLeft,xRight,yTop,yBottom,type):
    if not self.__staveFound:
      raise Exception("Cannot define a region for stave that has not been found.")
    
    if len(list(filter(lambda x : x > 1.0 or x < 0.0, [xLeft,xRight,yTop,yBottom]))) > 0:
      raise Exception("Regions are defined by relative coordinates w.r.t. the stave. The coordinates must be between 0.0 and 1.0.")
    
    if not(xLeft < xRight and yTop < yBottom):
      raise Exception("The coordiates for the region are invalid.")
    
    regionXLeft = self.__xLeft + xLeft*self.__length
    regionYTop = self.__yTop + yTop*self.__width
    regionXRight = self.__xLeft + xRight*self.__length
    regionYBottom = self.__yTop + yBottom*self.__width
    
    newRegion = Region(self.__globalImg,regionXLeft,regionXRight,regionYTop,regionYBottom)
    
    if type in self.__regions:
      self.__regions[type].append(newRegion)
    else:
	  self.__regions[type] = []
	  self.__regions[type].append(newRegion)
	
  def Echo(self):
    if self.__staveFound:
      print(str([self.__xLeft,self.__xRight,self.__yTop,self.__yBottom]))
    else:
      raise("Stave has not been found yet.")
    
   
  def PrintRegions(self,regionType):
    for region in self.__regions[regionType]:
	    region.Echo()
   
  def Show(self):
    plt.imshow(self.__globalImg)
    plt.show()
    return

  def SaveImage(self,path):
    plt.imshow(self.__globalImg)
    plt.savefig(path)
    return
	
  def DrawRegions(self,imgToBeImprinted,regionType):
    for region in self.__regions[regionType]:
      position = region.getPosition()
      cv2.rectangle(imgToBeImprinted,(int(position[0]),int(position[2])),(int(position[1]),int(position[3])),(0,0,0),thickness=1)
    return
    
  def DrawEdges(self,img):
    if self.__staveFound:
      cv2.rectangle(img,(self.__xLeft,self.__yTop),(self.__xRight,self.__yBottom),(0,0,0),thickness=1)
    else:
      raise Exception("Stave was not found. Cannot use DrawEdges.")
    return

  def getImage(self):
    deepCopy = np.copy(self.__globalImg)
    return deepCopy
  
  def getTemperatures(self,regionType):
    temperatures = []
    for region in self.__regions[regionType]:
      temperatures.append(region.getAverageTemperature())
    return temperatures
	
  def getImpedances(self,regionType):
    regionTemp = self.getTemperatures(regionType)
    #temp. of the cooling liquid
    liquidTemperature = self.__temperatureProfile
    
    #scale scale it up
    liquidTemperature = map(lambda x:x*((self.__Tout-self.__Tin)/self.__temperatureProfile[-1]), self.__temperatureProfile)

    #shift it to match Tin
    liquidTemperature = map(lambda x:x+(self.__Tin-liquidTemperature[0]), liquidTemperature)
    
    flowRateKgPerSec = self.__FR/60
    
    impedances = []

    for i in range(0,28):
      #divide by two to get the heat only for one part
      heat = (liquidTemperature[i] - liquidTemperature[i+1])*self.__heatCapacity*0.5*flowRateKgPerSec
      averageTempDiff = (liquidTemperature[i]+liquidTemperature[i+1])/2 - regionTemp[i]
      impedances.append(averageTempDiff/heat)
    
    return impedances
  
class Region:
  def __init__(self,globalImg,xLeft,xRight,yTop,yBottom):
    self.__xLeft = xLeft
    self.__xRight = xRight
    self.__yTop = yTop
    self.__yBottom = yBottom
    self.__img = globalImg[int(yTop):int(yBottom),int(xLeft):int(xRight)]
    self.__averageTemperature = np.mean(self.__img)
	
  def getAverageTemperature(self):
    return self.__averageTemperature
	
  def getPosition(self):
    return [self.__xLeft,self.__xRight,self.__yTop,self.__yBottom]
  
  def Echo(self):
    print(str([self.__xLeft,self.__xRight,self.__yTop,self.__yBottom]))
    


