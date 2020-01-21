#!/usr/bin/env python

'''
configFinder.py

Author: William Heidorn, Iowa State University
About: This program takes a thermal image of an ATLAS Itk Stave Support, and
  finds the stave and produces 4 points that contain the stave

Requires: pyROOT, Python 2.7, OpenCV

'''

import numpy as np
import cv2
import ROOT


def FindPoints(strImageFile,strOutputFile,outdir,bol14ModCore = False,xPixels = 640,yPixels = 480,fltxPercentCutL=0.05,fltxPercentCutR=0.023,fltyPercentCut=0.20):
  """
  This function takes an input root stave image and finds all of the appropriate
  locations on the stave and creates a config file.
  """

  try:
    #Load in the file
    imageFile = ROOT.TFile(strImageFile,"read")
  except:
    print("Failed to Load ImageFile")
    return
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



  # Make The Canny Image
  v = np.median(image)

  sigma = 0.33
  lower = int(max(0,(1-sigma)*v))
  upper = int(min(255,(1+sigma)*v))

  laplacian = cv2.Canny(np.uint8(image),lower,upper)
  image2 = laplacian

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


  # Find the Four Corners of the Pipe Area

  c2.lines = []  #The lines stored in the canvas to show the buildup
  lineData = []  #This will be the four points
  HorData  = []  #The average y value of each horizontal line
  VertData = []  #The average x value of each vertical line
  ShortHorData = []

  orighist.Draw("colz")
 
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
      lineObj = ROOT.TLine(x1,y1,x2,y2)
      lineObj.SetLineColor(3)
      lineObj.SetLineWidth(3)
      c2.lines += [lineObj]
      lineObj.Draw()

    HorData = np.sort(HorData) 
    CentSep = np.amax(abs(HorData-240))
    while CentSep > 50:
      lenHor = np.size(HorData)
      if abs(HorData[0]-240)> 50:
        HorData = np.delete(HorData,0)
      else:
        HorData = np.delete(HorData,lenHor-1)
      CentSep = np.amax(abs(HorData-240))

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
    StaveLength = 580 #~14 module stave core length (This is approximated from a Yale thermal image from Rec-000110)
    cutPercent= 0.02
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
        lineObj = ROOT.TLine(x1,y1,x2,y2)
        lineObj.SetLineColor(2)
        lineObj.SetLineWidth(4)
        c2.lines += [lineObj]
        lineObj.Draw()

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

  x0=lineData[3]
  x1=lineData[2]
  y0=lineData[1]
  y1=lineData[0]
  Dx = x1-x0 
  Dy = y1-y0
 
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
"""
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
