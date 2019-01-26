#!/usr/bin/python

"""
@run:
  ./Vignetting.py [files]

  Each of the files needs to be one of a metal block that is uniformly cooled
  and then shifted across the field of view in the long direction. There should
  be 15 data files which are in the attached folder

@brief:

@notes:

@email: wheidorn@iastate.edu
"""

import sys
import os
import ROOT
import numpy as np

def ReadInFrame(filename,outdir,stripnumber,stripdata,nc = 40, nx = 640, ny = 480):
  """
  This reads in an input root file and spits out a strip of data from it and
  writes it to the strip image
  """
  tfile = ROOT.TFile(filename,"read")
  atree = tfile.Get("atree")
  xpos = np.zeros(1,dtype = int)
  ypos = np.zeros(1,dtype = int)
  temp = np.zeros(1,dtype = float)

  atree.SetBranchAddress( "xpos", xpos)
  atree.SetBranchAddress( "ypos", ypos)
  atree.SetBranchAddress( "temperature", temp)
  nentries = atree.GetEntries()

  tempdata = [[0. for i in range(ny)] for j in range(nx)]
  for i in range(nentries):
    atree.GetEntry( i )
    tempdata[xpos[0]][ypos[0]] = temp[0]    

  #Print out the whole first strip's plot
  if stripnumber == 0:
    printPlot(tempdata,outdir,"FirstPlot")

  #Fill the old data into the output
  oldstripdata = stripdata
  outputdata = [[-999. for i in range(ny)] for j in range(nx)]
  for x in range(nx):
    for y in range(ny):
      outputdata[x][y] = oldstripdata[x][y]

  #Fill the strip into the output
  xf = nx - 20 - nc* stripnumber
  x0 = xf - nc  
  for x in range(x0,xf):
    for y in range(ny):
      outputdata[x][y] = tempdata[x][y]

  return outputdata
#------------------------------------------------------------------------------
def printPlot (dataArray,outdir,name = "StripTemps",nx = 640,ny = 480):
  """
  This takes in an array, fills a histogram and then plots it
  """
  ROOT.gStyle.SetOptStat(0)
  Hist = ROOT.TH2F('StripsHist',name+';Xpixels;Ypixels;Temperature',nx,0,nx,ny,0,ny)
  for x in range(nx):
    for y in range(ny):
      Hist.SetBinContent(x+1,y+1,dataArray[x][y])

  Canvas = ROOT.TCanvas("c1")
  Canvas.cd()
  Hist.Draw("colz")
  if name == "StripTemps" or name == "FirstPlot":
    Hist.GetZaxis().SetRangeUser(-40,20)
  else:
    Hist.GetZaxis().SetRangeUser(-45,-35)
  Canvas.Print(outdir+"/"+name+"Hist.png")

#------------------------------------------------------------------------------
def findXVig (dataArray, outdir,nstrips = 15,nc =40 ,nx =640,ny =480):
  """
  This takes in the total array and outputs 15 numbers which is the vignetting across the x direction
  """
  #The region of interest
  y0 = 180
  y1 = 220
  y2 = 260
  yf = 300
  

  npix = ((y1-y0)+(yf-y2))*nc

  xgraph = ROOT.TGraphErrors(15)
 
  outputData = [0 for x in range(nstrips)]
  outputError = [0 for x in range(nstrips)]

  #Find average in area
  for strip in range(nstrips):
    x0 = 20 + nc*strip
    xf = x0 + nc
    count = 0
    dataset = [0 for i in range(npix)]
    for x in range(x0,xf):
      for y in range(y0,y1): 
        dataset[count] = dataArray[x][y]
        count += 1
      for y in range(y2,yf):
        dataset[count] = dataArray[x][y]
        count += 1
    avg = np.average(dataset)
    std = np.std(dataset)
    outputData[strip]= avg
    outputError[strip] = std
  #print dataset
  #return
  for strip in range(nstrips):
    xgraph.SetPoint(strip,strip+1,outputData[strip] - outputData[7])
    xgraph.SetPointError(strip,0,(outputError[strip]+outputError[7])**0.5)
  canvas0 = ROOT.TCanvas()
  xgraph.SetTitle("Vignetting Across X; Strip Number [40pixel avg]; Average Vignetting [#circC]")
  xgraph.Draw("AP")
  xgraph.SetMarkerStyle(21)
  xgraph.SetMarkerColor(2)
  canvas0.SetGrid() 

  canvas0.Print(outdir+"/XdirectionVignetting.png") 


#------------------------------------------------------------------------------
def filterScrews (dataArray, outdir,nstrips = 15,nc =40 ,nx =640,ny =480):
  """
  This intakes an array, looks at each row and removes outlier points between them
  """
  #The region of interest
  y0 = 180
  yf = 300
  stdmult = 0.5

  npix = (300 - 180)*nc

  outputData = [[-999. for y in range(ny)] for x in range(nx)]

  #Find average in area
  for strip in range(nstrips):
    x0 = 20 + nc*strip
    xf = x0 + nc
    count = 0
    dataset = [0 for i in range(npix)]
    for x in range(x0,xf):
      for y in range(y0,yf):
        dataset[count] = dataArray[x][y]
        count += 1
    avg = np.average(dataset)
    std = np.std(dataset)
    minT = avg - std*stdmult
    maxT = avg + std*stdmult

    for x in range (x0,xf):
      for y in range(y0,yf):
        Temp = dataArray[x][y]
        if Temp > minT and Temp < maxT: 
          outputData[x][y] = dataArray[x][y]
        else:
          outputData[x][y] = avg  

  printPlot(outputData,outdir,"Filtered")
  return outputData

#------------------------------------------------------------------------------
def printVignetting (dataArray,outdir,nstrips = 15,nc = 40,nx = 640,ny = 480):
  """
  This takes in an array, and condenses it to an average along the y direction
  """
  AvgedData = [[0. for y in range(ny)] for x in range(nstrips)]
  AvgHist = ROOT.TH2F('AvgHist','Avg Temps in x in Each Strip;Xpixels;Ypixels;Temperature',nstrips,20,nx-20,ny,0,ny)
  
  y0 = 180
  yf = 300
  ny2 = yf - y0

  #This is the blank vignetting temperature map
  VHist = ROOT.TH2F('VHist','Vignetting Temps;Xpixels;Ypixels;Temperature Diff',nstrips,20,nx-20,ny2,y0,yf)

  #Fill the map from the strips
  for strip in range(nstrips):
    x0 = 20 + nc*strip
    xf = x0 + nc
    for y in range(ny):
      for x in range(x0,xf):
        #Average over the x direction in the strips
        AvgedData[strip][y] += dataArray[x][y]/nc
  for y in range(ny):
    for strip in range(nstrips):
      AvgHist.SetBinContent(strip+1,y+1,AvgedData[strip][y])

  for y in range(ny2):
    for strip in range(nstrips):
      VHist.SetBinContent(strip+1,y+1,AvgedData[strip][y+y0]-AvgedData[7][y+y0])

  maxVig = 0
  avgVig = 0
  #Create an abs value vignetting plot
  VHistAbs = VHist.Clone()
  for y in range(ny2):
    for strip in range(nstrips):
      VHistAbs.SetBinContent(strip+1,y+1,abs(AvgedData[strip][y+y0]-AvgedData[7][y+y0]))
      maxVig = max(maxVig,abs(AvgedData[strip][y+y0]-AvgedData[7][y+y0]))
      avgVig += abs(AvgedData[strip][y+y0]-AvgedData[7][y+y0])/(ny2*nstrips)

  print("Max Vignetting: "+str(maxVig))
  print("Avg Vignetting: "+str(avgVig))

  #Create a condensed vignetting plot
  nYstrips=11
  FinalV = [[0. for y in range(nYstrips)] for x in range(nstrips)]
  VHistF = ROOT.TH2F('VHistF','Vignetting ;;;Temperature Diff',nstrips,20,nx-20,nYstrips,20,ny-20)
  for ystrip in range(nYstrips):
    for x in range(nstrips):
      y0 =20 + nc*ystrip
      yf = y0+nc
      for y in range(y0,yf):
        FinalV[x][ystrip] += AvgedData[x][y]/40

  for ystrip in range(nYstrips):
    for x in range(nstrips):
      VHistF.SetBinContent(x+1,ystrip+1,FinalV[x][ystrip]-FinalV[7][5])

  #Print out results
  Canvas = ROOT.TCanvas("c2")
  Canvas.cd()
  ROOT.gStyle.SetOptStat(0)
  AvgHist.Draw('colz')
  AvgHist.GetZaxis().SetRangeUser(-45,-35)
  Canvas.Print(outdir+'/AvgHist.png')   
  ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)
  VHist.Draw('colz')
  VHist.SetContour(100)
  VHist.GetZaxis().SetRangeUser(-5,5)
  Canvas.Print(outdir+'/VHist.png')
  VHistAbs.Draw('colz')
  VHistAbs.SetContour(100)
  VHistAbs.GetZaxis().SetRangeUser(-0.00001,5)
  Canvas.Print(outdir+'/VHistAbs.png')
  VHistF.Draw('cont4z')
  VHistF.GetYaxis().SetRangeUser(y0,yf)
  Canvas.Print(outdir+'/VhistF.png')

#------------------------------------------------------------------------------
def main ():
  """
  The main routine
  """
  if sys.version_info[0] >=3:
    print ("ERROR: PyROOT only works for python 2.x")
    raise Exception(" Python Version too high.")

  #Load the filenames
  nargv = len(sys.argv)
  inputfiles = []
  for i in range(1,nargv):
    inputfiles = np.append(inputfiles,sys.argv[i])

  nfiles = len(inputfiles)

  nX = 640
  nY = 480
  LoadedDataStrips = [[-999. for i in range(nY)] for j in range(nX)]

  outdir = inputfiles[0]
  outdir = outdir.split('/')[0]

  for i in range(nfiles):
    LoadedDataStrips = ReadInFrame(inputfiles[i],outdir,i,LoadedDataStrips)

  printPlot(LoadedDataStrips,outdir)
  findXVig (LoadedDataStrips, outdir)
  LoadedDataStrips = filterScrews(LoadedDataStrips,outdir)
  printVignetting(LoadedDataStrips,outdir)

if __name__ == "__main__":
  main()
