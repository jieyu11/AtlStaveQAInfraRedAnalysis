#!/usr/bin/env python

"""
stave.py

Author: Lubos Vozdecky (Queen Mary, University of London)
About: This code implements the Stave and Region class. That are used in the impedanceFromCSV.py script. It uses a simple algorithm
for finding the edges of the staves.
"""

import logging
import os
import cv2
import configparser
import numpy as np
from matplotlib import pyplot as plt

class Stave:
  
  def __init__(self,globalImg, configFile):
    self.__globalImg = globalImg
    self.__staveFound = False
    self.__xLeft = 0
    self.__xRight = 0
    self.__yTop = 0
    self.__yBottom = 0
    self.__regions = {} # dictionary
    self.__staveRatio = 11.957 #length:width ratio for sanity check
    self.__staveRatioTolerance = 0.05 # 5% tolerance on the stave length/width ratio
    self.__lineThickness = 2 #thickness of the line that is used for drawing the regions
    

    #importing data from the config file
    if not os.path.isfile(configFile):
      raise Exception("Could not find config file '" + configFile + "'")
	
    logging.debug("Importing variables from config file " + configFile + ":")
    config = configparser.ConfigParser()
    config.read(configFile)

    self.__Tin = float(config["Default"]["temp_in"])
    self.__Tout = float(config["Default"]["temp_out"])
    self.__heatCapacity = float(config["Default"]["c_liquid"])
    self.__FR = float(config["Default"]["flow_rate"])
    self.__temperatureProfile = [float(x) for x in config["Default"]["temperatureProfile"].split(",")]
    
    #print the imported variables
    for variable in config.items("Default"):
      logging.debug(variable[0] + " = " + variable[1])
    
  def FindStaveWithin(self, relXMin, relXMax, relYMin, relYMax):
    if relXMin > relXMax or relYMin > relYMax:
      logging.error("Minimal values are larger than the maximum.")
      raise Exception("FindStaveWithin(): minimal values are larger than the maximum.")
    
    if self.__staveFound:
      logging.error("Calling FindStaveWithin on a stave that was already found.")
      raise Exception("Stave was already found. Cannot call FindStaveWithin again.")
    
    logging.debug("Finding stave algorithm:")
    
    imageHeight = self.__globalImg.shape[0]
    imageWidth = self.__globalImg.shape[1]
    
    xMin = int(relXMin * imageWidth)
    xMax = int(relXMax * imageWidth)
    yMin = int(relYMin * imageHeight)
    yMax = int(relYMax * imageHeight)
    
    imgOfInterest = self.__globalImg[yMin:yMax,xMin:xMax]
    logging.debug("Looking for the stave within: " + str([yMin,yMax,xMin,xMax]))
    
    gradient = np.gradient(imgOfInterest,axis=0)
    
    projection = np.sum(gradient, axis=1)
    
    upperEdge = np.argmax(projection)
    logging.debug("upperEdge = " + str(upperEdge))
    lowerEdge = np.argmin(projection)
    logging.debug("lowerEdge = " + str(lowerEdge))
    
    strip = imgOfInterest[upperEdge:lowerEdge,:]
    
    stripGradient = np.gradient(strip,axis=1)
    
    stripProjection = np.sum(stripGradient, axis=0)
    
    leftEdge = np.argmax(stripProjection)
    logging.debug("leftEdge = " + str(leftEdge))
    rightEdge = np.argmin(stripProjection)
    logging.debug("rightEdge = " + str(rightEdge))
    
    ratioMeasured = abs(1.0*(rightEdge-leftEdge)/(lowerEdge-upperEdge))
    logging.debug("ratioMeasured = " +str(ratioMeasured))
    
    
    if abs(self.__staveRatio - ratioMeasured)/self.__staveRatio > self.__staveRatioTolerance:
      logging.error("Ratio test FAILED")
      raise Exception("Ratio of the found stave is not within the limits.")
    
    logging.debug("Ratio test PASSED")
    
    self.__xLeft = leftEdge + xMin
    self.__xRight = rightEdge + xMin
    self.__yTop = upperEdge + yMin
    self.__yBottom = lowerEdge +yMin
    
    self.__width = self.__yBottom - self.__yTop
    self.__length = self.__xRight - self.__xLeft

    self.__staveFound = True
    
    
  def ScaleImage(self, scale):
    logging.debug("Scaling image: scale = " + str(scale))
    width = int(self.__globalImg.shape[1] * scale)
    height = int(self.__globalImg.shape[0] * scale)
    dim = (width, height)
    resized = cv2.resize(self.__globalImg, dim, interpolation = cv2.INTER_LINEAR)
    self.__globalImg = resized
    
    self.__lineThickness = self.__lineThickness * scale
  
  def setTemperatureProfile(self,newTemperatureProfile):
    logging.debug("Re-setting the temperature profile:")
    self.__temperatureProfile = newTemperatureProfile
    logging.debug("temperatureProfile = " + str(self.__temperatureProfile))

  def AddRegion(self,xLeft,xRight,yTop,yBottom,type):
    if not self.__staveFound:
      logging.error("Defining region for a stave that has not been found.")
      raise Exception("Cannot define a region for stave that has not been found.")
    
    if len(list(filter(lambda x : x > 1.0 or x < 0.0, [xLeft,xRight,yTop,yBottom]))) > 0:
      logging.error("Invalid coordinates for a region. [xLeft,xRight,yTop,yBottom] = " + str([xLeft,xRight,yTop,yBottom]))
      raise Exception("Regions are defined by relative coordinates w.r.t. the stave. The coordinates must be between 0.0 and 1.0.")
    
    if not(xLeft < xRight and yTop < yBottom):
      logging.error("Invalid coordinates for a region. [xLeft,xRight,yTop,yBottom] = " + str([xLeft,xRight,yTop,yBottom]))
      raise Exception("The coordiates for the region are invalid.")
    
    logging.debug("Adding a new region of type '" + str(type) + "': [xLeft,xRight,yTop,yBottom] = " + str([xLeft,xRight,yTop,yBottom]))
    
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
	
  def AddUBendRegion(self, rxLeft, rxRight, ryTop, ryBottom, rradius, rlength, type, bend="downwards"):
    if not self.__staveFound:
      logging.error("Defining region for a stave that has not been found.")
      raise Exception("Cannot define a region for stave that has not been found.")
    
    if len(list(filter(lambda x : x > 1.0 or x < 0.0, [rxLeft, rxRight, ryTop, ryBottom]))) > 0:
      logging.error("Invalid coordinates for a region. [xLeft,xRight,yTop,yBottom] = " + str([rxLeft, rxRight, ryTop, ryBottom]))
      raise Exception("Regions are defined by relative coordinates w.r.t. the stave. The coordinates must be between 0.0 and 1.0.")
    
    if not(rxLeft < rxRight and ryTop < ryBottom):
      logging.error("Invalid coordinates for a region. [xLeft,xRight,yTop,yBottom] = " + str([rxLeft, rxRight, ryTop, ryBottom]))
      raise Exception("The coordiates for the region are invalid.")
    
    if not(bend=="downwards" or bend=="upwards"):
      logging.error("Invalid  argument bend = " + bend)
      raise Exception("Invalid  argument bend.")
    
    logging.debug("Adding a new region of type '" + str(type) + "': [xLeft,xRight,yTop,yBottom] = " + str([rxLeft, rxRight, ryTop, ryBottom]))
    
    #coordinates of the first "starting" rectangle
    xLeft1 = int(self.__xLeft + rxLeft*self.__length)
    yTop1 = int(self.__yTop + ryTop*self.__width)
    xRight1 = int(self.__xLeft + rxRight*self.__length)
    yBottom1 = int(self.__yTop + ryBottom*self.__width)
    
    
    radius = int(rradius*self.__width)
    
    thickness = int(yBottom1-yTop1)
    centreX = int(xRight1)
    
    if bend=="upwards":
      centreY = int(0.5*(yTop1+yBottom1)-rradius*self.__width)
    else:
      centreY = int(0.5*(yTop1+yBottom1)+rradius*self.__width)
    
    xLeft2 = int(xRight1 + radius - thickness/2)
    xRight2 = int(xLeft2 + thickness)
    
    if bend=="upwards":
      yTop2 = int(yBottom1 - radius - thickness/2)
      yBottom2 = int(yTop2 - rlength*self.__width)
    else:
      yTop2 = int(yBottom1 + radius - thickness/2)
      yBottom2 = int(yTop2 + rlength*self.__width)
    
    
    if bend=="upwards":
      start_angle, stop_angle = 0, 90
    else:
      start_angle, stop_angle = 270, 360
    
    shape = np.shape(self.__globalImg)
    
    regions_image = np.zeros([shape[0],shape[1]])

    regions_image = cv2.rectangle(regions_image, (xLeft1,yTop1), (xRight1,yBottom1), 1, -1) 

    regions_image=cv2.ellipse(regions_image, (centreX,centreY), (radius,radius), 0, start_angle, stop_angle, 1, thickness)

    regions_image = cv2.rectangle(regions_image, (xLeft2,yTop2), (xRight2,yBottom2), 1, -1) 
    
    newRegion = GeneralRegion(self.__globalImg, regions_image)
    
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
      region.DrawRegion(imgToBeImprinted,self.__lineThickness)
    return
    
  def DrawEdges(self,img):
    if self.__staveFound:
      cv2.rectangle(img,(self.__xLeft,self.__yTop),(self.__xRight,self.__yBottom),(0,0,0),thickness=self.__lineThickness)
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
    logging.debug("Calculating impedances of region type '" + str(regionType) + "'")
    regionTemp = self.getTemperatures(regionType)
    #temp. of the cooling liquid
    liquidTemperature = self.__temperatureProfile
    
    logging.debug("Scaling the temperature profile accoreding to [Tout,Tin] = " + str([self.__Tout,self.__Tin]))
    
    #scale it up
    liquidTemperature = [x*(self.__Tout-self.__Tin)/self.__temperatureProfile[-1] for x in self.__temperatureProfile]

    #shift it to match Tin
    liquidTemperature = [x+(self.__Tin-liquidTemperature[0]) for x in liquidTemperature]
    
    flowRateKgPerSec = self.__FR/60
    
    logging.debug("liquidTemperature after scaling = " + str(liquidTemperature))
    logging.debug("flowRateKgPerSec = " + str(flowRateKgPerSec))
    logging.debug("heatCapacity = " + str(self.__heatCapacity))
    
    impedances = []
    
    if not len(regionTemp)+1==len(liquidTemperature):
      logging.error("The number of temperature profile data points does not match the number of regions.")
      raise Exception("The number of temperature profile data points does not match the number of regions.")

    for i in range(0,len(regionTemp)):
      #divide by two to get the heat only for one part
      heat = abs(liquidTemperature[i] - liquidTemperature[i+1])*self.__heatCapacity*0.5*flowRateKgPerSec
      averageTempDiff = abs((liquidTemperature[i]+liquidTemperature[i+1])/2 - regionTemp[i])
      impedances.append(averageTempDiff/heat)
    
    logging.debug("impedances = " +str(impedances))
    
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
    
  def DrawRegion(self, imgToBeImprinted, thickness):
    cv2.rectangle(imgToBeImprinted,(int(self.__xLeft),int(self.__yTop)),(int(self.__xRight),int(self.__yBottom)),(0,0,0),thickness=thickness)
    

#added for implementing the U-bend regions
class GeneralRegion:
  def __init__(self, globalImg, regionImg):
    self.__globalImg = globalImg
    self.__regionImg = regionImg
    
    selected_regions = globalImg*regionImg
    self.__selected_regions = selected_regions
    totalSum = np.sum(selected_regions)
    numberNonzero = len(np.nonzero(selected_regions)[0])

    self.__averageTemperature = totalSum/numberNonzero
    
  def getAverageTemperature(self):
    return self.__averageTemperature
    
  def getPosition(self):
    return [np.min(np.nonzero(self.__selected_regions)[1]),np.max(np.nonzero(self.__selected_regions)[1]), np.min(np.nonzero(self.__selected_regions)[0]),np.max(np.nonzero(self.__selected_regions)[0])]
    
  def DrawRegion(self, imgToBeImprinted, thickness):
    #update the value, not the reference -> use [:]
    imgToBeImprinted[:] = imgToBeImprinted + 100*self.__regionImg