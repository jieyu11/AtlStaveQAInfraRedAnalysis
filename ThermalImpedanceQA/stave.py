#!/usr/bin/env python

import logging
import os
import cv2
import configparser
import numpy as np
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
    
  def FindStaveWithin(self, xMin, xMax, yMin, yMax):
    if xMin > xMax or yMin > yMax:
      raise Exception("FindStaveWithin(): minimal values are larger than the maximum.")
    
    if self.__staveFound:
      raise Exception("Stave was already found. Cannot call FindStaveWithin again.")
    
    
    imgOfInterest = self.__globalImg[yMin:yMax,xMin:xMax]
    
    xPixels = imgOfInterest.shape[0]
    yPixels = imgOfInterest.shape[1]
    
    # Make The Canny Image
    v = np.median(imgOfInterest)
    sigma = 0.33
    lower = int(max(0,(1-sigma)*v))
    upper = int(min(255,(1+sigma)*v))

    cannyImage = cv2.Canny(np.uint8(imgOfInterest),lower,upper)
  
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
      findLongLines = cv2.HoughLinesP(cannyImage,rho = 1,theta = 1*np.pi/1000,threshold = 100,minLineLength = 100,maxLineGap = 10)
      findLongLines = findLongLines.flatten()
      nlines = len(findLongLines)/4
      findLongLines = np.split(findLongLines,nlines)
      
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
      
      upperEdge = np.max(np.extract(HorData<=aproxYpos,HorData)[0])      
      lowerEdge = np.min(np.extract(HorData>aproxYpos,HorData)[0])

      lineData += [lowerEdge]
      lineData += [upperEdge]

    except:
      logging.error("Failed to Find Long Horizontal Lines")
      return
    
    
    
    StaveLength = 630 #~14 module stave for the QMUL IR set-up
    cutPercent = 0.20 # QMUL
    
    #Find Short Vert Lines
    try:
      findShortLines = cv2.HoughLinesP(cannyImage,rho = 1,theta = 1*np.pi/10000,threshold = 20,minLineLength = 10, maxLineGap = 5)

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

  def AddRegion(self,xLeft,xRight,yTop,yBottom,type):
    if not self.__staveFound:
	  raise Exception("Cannot define a region for stave that has not been found.")
	
    if type in self.__regions:
      self.__regions[type].append(Region(self.__staveImg,xLeft,xRight,yTop,yBottom))
    else:
	  self.__regions[type] = []
	  self.__regions[type].append(Region(self.__staveImg,xLeft,xRight,yTop,yBottom))
	
  def Echo(self):
    if self.__staveFound:
      print(str([self.__xLeft,self.__xRight,self.__yTop,self.__yBottom]))
    else:
      raise("Stave has not been found yet.")
    
   
  def PrintRegions(self):
    for regionType in self.__regions:
      print(regionType+": ")
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
      cv2.rectangle(imgToBeImprinted,(position[2],position[3]),(position[0],position[1]),(0,0,0),thickness=1)
    return
    
  def DrawEdges(self,img):
    if self.__staveFound:
      cv2.rectangle(img,(self.__xLeft,self.__yTop),(self.__xRight,self.__yBottom),(0,0,0),thickness=1)
    else:
      raise("Stave was not found. Cannot use DrawEdges.")
    return
	   
class Region:
  def __init__(self,staveImg,xLeft,xRight,yTop,yBottom):
    self.__xLeft = xLeft
    self.__xRight = xRight
    self.__yTop = yTop
    self.__yBottom = yBottom
    self.__img = staveImg[xLeft:xRight,yTop:yBottom]
    self.__averageTemperature = np.mean(self.__img)
	
  def getAverageTemperature():
    return self.__averageTemperature
	
  def getPosition():
    return [self.__xLeft,self.__xRight,self.__yTop,self.__yBottom]
  
  def Echo(self):
    print(str([self.__xLeft,self.__xRight,self.__yTop,self.__yBottom]))
    


