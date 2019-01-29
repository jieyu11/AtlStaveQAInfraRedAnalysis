#!/usr/bin/python
"""
@run
  ./SliceFinder.py Path/To/frame_average.root

  parameters

@brief
  This code takes an input temperature file in either csv form or ROOT form.
  It then converts it using methods described by BNL into an Impedance Map
  across the stave. A Basic outline of the conversion is given below.

  0.0 Read in image file, attempt to read in a config_file that has the 
    stave pixel cut. It then cuts the image to the size of the stave.

  0.1 Finds the pixel to m conversion

  1.0 The image is smoothed to remove shiny spots caused by the metal spots
    on the stave.
@notes

  This program takes a frame_average.root file and converts it into average 
  temperature strips across each of the module locations. It loads in a csv
  file named FEATemps.csv

@reference

@email
  wheidorn@iastate.edu
"""
import sys
import os
import ROOT
import numpy as np
from scipy.signal import argrelextrema

#To get the config in the event it doesn't exist... This currently doesn't work
sys.path.append("../configFinder.py")
#import configFinder as cf

#------------------------------------------------------------------------------
#LOADING IN THINGS
#------------------------------------------------------------------------------
def ReadInCSV(strInputFile,strName,strOutdir):
  """
  Takes in a CSV file and then outputs the data as a nested list
  """
  csvData = open(strInputFile,'r')
  lineCount=0
  dataCount=0
  for line in csvData:
    lineCount+=1
    ldata = line.split(',')
    dataCount = len(ldata)

  ny = lineCount
  nx = dataCount
  lstTempData = [[0 for y in range(ny)] for x in range(nx)]

  csvData.close()
  csvData = open(strInputFile,'r')

  count = 0
  for line in csvData:
    strline = line
    strline = strline.strip('\r\n')
    strline = strline.strip(' ')
    fltline = strline.split(',') 
    fltline = [float(i) for i in fltline]
    for x in range(nx):
      lstTempData[x][count] = fltline[x]
    count+=1
  csvData.close()
 
  CutHist = ROOT.TH2F("h2",strName+" Loaded Image;xPixel;yPixel",nx,0,nx,ny,0,ny) 
  for x in range(nx):
    for y in range(ny):
      CutHist.SetBinContent(x+1,y+1,lstTempData[x][y])

  c0 = ROOT.TCanvas("c0")
  c0.cd()
  CutHist.Draw("COLZ")
  c0.Print(strOutdir+strName+"_LoadedImage.png")
  return lstTempData

#------------------------------------------------------------------------------
def CutCSV(lstTempData,strName,strOutdir,X0=51,X1=444,Y0=8,Y1=198):
  """
  This cuts an input CSV file
  """
  nyCut = Y1- Y0 
  nxCut = X1- X0

  lstCutTempData = [[0 for y in range(nyCut)] for x in range(nxCut)]
  for x in range(nxCut):
    for y in range(nyCut):
      lstCutTempData[x][y] = lstTempData[x+X0][y+Y0]

  CutHist = ROOT.TH2F("h2",strName+" Cut Image;xPixel;yPixel",nxCut,0,nxCut,nyCut,0,nyCut) 
  for x in range(nxCut):
    for y in range(nyCut):
      CutHist.SetBinContent(x+1,y+1,lstCutTempData[x][y])

  c0 = ROOT.TCanvas("c0")
  c0.cd()
  CutHist.Draw("COLZ")
  c0.Print(strOutdir+strName+"_Cut_Image.png")
  
  return [lstCutTempData,nxCut,nyCut,58-Y0,136-Y0]

#------------------------------------------------------------------------------
def LoadSlices(inputcsv,outdir):
  """
  reads in a csv with strip data and outputs an array
  """
  csvData = open(inputcsv,'r')
  lineCount = 0
  dataCount = 0
  for line in csvData:
    lineCount +=1
    ldata = line.split(',')
    dataCount = max(len(ldata),dataCount)

  csvData.close()
  nx = dataCount
  ny = lineCount

  Data = [[-999 for y in range(ny)] for x in range(nx)]

  count = 0
  MaxT = -55
  MinT = 55

  csvData = open(inputcsv,'r')
  for line in csvData:
    strline = line 
    strline = strline.strip('\r\n')
    strline = strline.strip(' ') 
    strline = strline.split(',')
    fltline = [0 for i in range(len(strline))]
    #print strline
    for i in range(len(strline)):
      if '-' in strline[i]:
        #print type(strline[i])
        num = float(strline[i].split('-')[-1]) 
        num = -num
      elif '' in strline[i]:
        num = -999
      else:
        num = float(strline[i])
      fltline[i] = num
    #print fltline 
    for x in range(nx):
      Data[x][count] = fltline[nx-x-1]
      if fltline[x]!= -999:
        MaxT = max(fltline[x],MaxT)
        MinT = min(fltline[x],MinT)
    count += 1
  csvData.close()

  csvHist = ROOT.TH2F("csvHist","Loaded Temp Lines; Module Number; yPosition [mm]",nx,0,nx,ny,0,ny)
  for x in range(nx):
    for y in range(ny): 
      csvHist.SetBinContent(x+1,y+1,Data[x][y])

  c0 = ROOT.TCanvas("c0","",8000,300)
  c0.cd()
  csvHist.Draw("COLZ")
  csvHist.GetZaxis().SetRangeUser(-40,0)
  c0.Print(outdir+"LoadedTempProfile.png")
  return [Data,nx,ny,.001,0.001]


#------------------------------------------------------------------------------
def ChopTempData(strInputFile,bolStaveSideL,X0,X1,Y0,Y1,strName,strOutdir):
  """
  Takes and loads the TTree from the input file and creates an array that has been cut to the size
  of the stave. It also will remove the End of Stave Card
  """
  File = ROOT.TFile(strInputFile,'r')
  Tree = File.Get('atree;1')
  nentries = Tree.GetEntries()
  ROOT.gStyle.SetOptStat(0)
  strName = strInputFile.split('/')[-2]
  print ("FILE NAME: "+strName)

  nxpixels = 640
  nypixels = 480
  Temp = [[0 for i in range(nypixels)] for j in range(nxpixels)]

  if bolStaveSideL <= 0:
    for y in range(nypixels):
      for x in range(nxpixels):
        Tree.GetEntry(x*nypixels + y)
        Temp[x][y] = Tree.temperature

  else: #Flipy along y direction!
    for y in range(nypixels):
      for x in range(nxpixels):
        Tree.GetEntry(x*nypixels + y)
        Temp[x][nypixels-1-y] = Tree.temperature

  nxcut = X1 - X0
  nycut = Y1 - Y0
  Tempcut = [[0 for i in range(nycut)] for j in range(nxcut)]
  for y in range(nycut):
    for x in range(nxcut):
      Tempcut[x][y] = Temp[x+X0][y+Y0]

  YEOScut = 5*nycut/7
  EOScut = [[0 for i in range(YEOScut)] for j in range(nxcut)]

  for y in range(YEOScut):
    for x in range(nxcut):
      EOScut[x][y] = Temp[x+X0][y+Y0]

  CutHist = ROOT.TH2F("h2",strName+" Cut Image;xPixel;yPixel;Temperature [#circC]",nxcut,0,nxcut,nycut,0,nycut) 
  for x in range(nxcut):
    for y in range(nycut):
      CutHist.Fill(x+1,y+1,Tempcut[x][y])

  c0 = ROOT.TCanvas("c0","",8000,300)
  c0.cd()
  CutHist.Draw("COLZ")
  CutHist.GetZaxis().SetRangeUser(-40,0)

  c0.Print(strOutdir+strName+"_CutImage.png")
  PipeB = 15
  PipeEOS = 35
  print nycut
  return [Tempcut,nxcut,nycut,PipeB,PipeEOS]
  #return [EOScut,nxcut,YEOScut,PipeB,PipeEOS]

#------------------------------------------------------------------------------
def RemoveShine(lstCutTempData,strName,strOutdir):
  """
  Smooths the areas of the stave that have shining spots
  """
  nx = lstCutTempData[1]
  ny = lstCutTempData[2]
  SmoothedData = lstCutTempData[0]

  SmoothBot = int(ny*0.15)
  SmoothTop = int(ny*0.17)
  #To Be figured out...
  for y in range(ny):
    if y > SmoothBot and y < ny-SmoothTop:
    #!!!if y > 40 and y < ny-30:
    #!!!if y > 7 and y < ny - 6:
    #!!!if y > 0:
      #print ('pixLine '+str(y)+" Not Smoothed")
      continue
    OrigLine = [0 for x in range(nx)]
    NewLine = [-999 for x in range(nx)]
    FinalLine = [-999 for x in range(nx)]
    
    #Read out the original line from the data set
    for x in range(nx):
      OrigLine[x] = SmoothedData[x][y] 

    #Find mean/stddev and remove things outside stdev
    mean = np.mean(OrigLine)
    stdDev = np.std(OrigLine)*0.5
    
    for x in range(nx):
      Temp = OrigLine[x]    
      if Temp > mean + stdDev:
        continue
      if Temp < mean - stdDev:
        continue
      NewLine[x] = Temp
      FinalLine[x] = Temp
   
    BolBadPoints = True
    while BolBadPoints == True:
      if -999 in NewLine:
        NewLine.remove(-999)
      else:
        BolBadPoints = False
    Nmean = np.mean(NewLine)

    #Put back the new mean value
    for x in range(nx):
      Temp = FinalLine[x]
      if Temp == -999:
        FinalLine[x] = Nmean

    for x in range(nx):
      SmoothedData[x][y] = FinalLine[x]
 
  CutHist = ROOT.TH2F("h2",strName+"_Smoothed Data;xPixel;yPixel",nx,0,nx,ny,0,ny) 
  for x in range(nx):
    for y in range(ny):
      CutHist.SetBinContent(x+1,y+1,SmoothedData[x][y])

  #Find Minima/Maxima of the plot
  EOSPipe = [-999 for i in range(nx)]
  BotPipe = [-999 for i in range(nx)]
  CentStave = [-999 for i in range(nx)]
  for x in range(nx):
    EOSYmax = -999
    BotYmax = -999
    for y in range(ny/2):
      if abs(SmoothedData[x][y]) > EOSYmax:
        EOSPipe[x] = y
        EOSYmax = abs(SmoothedData[x][y])
      if abs(SmoothedData[x][y+ny/2]) > BotYmax:
        BotPipe[x] = y + int(ny/2)
        BotYmax = abs(SmoothedData[x][y+ny/2])
    for x in range(nx):
      CentStave[x] = (BotPipe[x]+EOSPipe[x])/2

  c0 = ROOT.TCanvas("c0")
  c0.cd()
  CutHist.Draw("COLZ")
  c0.Print(strOutdir+strName+"_Smoothed.png")

  return [SmoothedData,nx,ny,EOSPipe,BotPipe,CentStave]

#------------------------------------------------------------------------------
def LoadStave( strInputFile, strName, strOutdir):
  """
  This takes an input file name and produces a plot of it's impedance and outputs
  an array of the impedance data
  """
  print('\nBEGIN TO LOAD STAVE FROM: '+strInputFile)

  bolDataFromCSV = False
  if '.csv' in strInputFile:
    print("Assuming Data Cut to Size Already")
    bolDataFromCSV = True

   #0 Read in the basic data and cut to stave body
  if bolDataFromCSV == False:
    #Load in the config file
    X0 = -999
    X1 = -999
    Y0 = -999
    Y1 = -999
    bolStaveSideL = -999

    ActualFileName = str(strInputFile.split('/')[-1])
    PossConfigLoc = strInputFile.replace(ActualFileName,'config_frame') 
    try:
      config_frame = open(PossConfigLoc,'r')
    except: 
      print("No config file found here... Creating a temporary one. For some reason the config Finder causes a segment violation...")
      Info = cf.FindPoints(strInputFile,'temp_config_frame')
      config_frame = open('temp_config_frame','r')
     
    for line in config_frame:
      try:
        Val = line.split(' ')
        Val = Val[-1]
        Val = Val.strip('\n')
        Val = int(Val)
      except: continue
      if 'SideL' in line:
        bolStaveSideL = Val
      if 'StavePixel' not in line: continue
      if 'X0' in line:
        X0 = Val
      elif 'X1' in line:
        X1 = Val
      elif 'Y0' in line:
        Y0 = Val
      elif 'Y1' in line:
        Y1 = Val
    print("FROM CONFIG FILE: X0 = {0} X1 = {1} Y0 = {2} Y1 = {3} SideL = {4} ".format(X0,X1,Y0,Y1,bolStaveSideL))  
 
    ChoppedTempData = ChopTempData(strInputFile,bolStaveSideL,X0,X1,Y0,Y1,strName,strOutdir)
    config_frame.close()
  else:    
    LoadedData = ReadInCSV(strInputFile,strName,strOutdir)
    ChoppedTempData = CutCSV(LoadedData,strName,strOutdir)

  #0.1 Get Pixel to m constants
  if bolDataFromCSV == True:
    fltPixToMeterX = 0.00328
    fltPixToMeterY = 0.000588
  else:
    fltPixToMeterX = 0.0023272727
    fltPixToMeterY = 0.0023272727
  print('LOADED FILE. SMOOTHING PROBLEM AREAS...')

  #1 Smooth out the shiny spots... Needs improvement
  SmoothedTempData = ChoppedTempData#RemoveShine(ChoppedTempData,strName,strOutdir)
 
  return [SmoothedTempData[0],SmoothedTempData[1],SmoothedTempData[2],fltPixToMeterX,fltPixToMeterY]

#------------------------------------------------------------------------------
def GetSlices(DataSet,intSlices,outdir,cutPercent =0.45):
  """
  This reads in a set of data and divides it into intSlices areas in the y dir
  it then takes the middle 80% and averages them to a single slice and then
  outputs a dataset with intSlices of data
  """
  Data = DataSet[0]
  nx = DataSet[1]
  ny = DataSet[2] 
  nxmod = int(nx/intSlices)
  nx0cut = int(nxmod*cutPercent)
  nxavg = int(nxmod - 2*nx0cut)

  print("Dividing nx = "+str(nx)+" into "+str(intSlices)+" averaged over "+str(nxavg))

  lineData = [[0. for y in range(ny)] for x in range(intSlices)]

  print("nx = "+str(nx)+" should = "+str(len(Data)))
  print("ny = "+str(ny)+" should = "+str(len(Data[0])))

  for i in range(intSlices):
    for y in range(ny):
      for x in range(nxavg):
        #print(str(x)+' '+str(y)+' '+str(i))
        #print str(x+nx0cut+nxmod*i) 
        lineData[i][y] += float(Data[x+nx0cut+nxmod*i][y])/nxavg

  avgHist = ROOT.TH2F("avgHist","Avg Temp Lines; Module Number; yPosition [m]",intSlices,0,intSlices,ny,0,ny*DataSet[4])
  for x in range(intSlices):
    for y in range(ny): 
      avgHist.SetBinContent(x+1,y+1,lineData[x][y])

  c0 = ROOT.TCanvas("c0","",8000,300)
  c0.cd()
  avgHist.Draw("COLZ")
  avgHist.GetZaxis().SetRangeUser(-40,0)
  c0.Print(outdir+"AvgTempProfile.png") 
  return [lineData,intSlices,ny,DataSet[3],DataSet[4]]
#------------------------------------------------------------------------------
def FindMinMaxAsy(Data,PixelCon):
  """
  From an input list and pixel conversion, this returns a list of information 
    [minima bottom,minima EOS,central maxima,asymmetry between minima]
  """
  npArray = np.array(Data)
  maxInd = argrelextrema(npArray,np.greater)[0]
  minInd = argrelextrema(npArray,np.less)[0]
  maxima = npArray[maxInd]
  minima = npArray[minInd]

  CMax = -999.
  CPCount = 0
  for i in range(len(maxInd)):
    Point = maxInd[i]*PixelCon
    if Point > .04 and Point < .06 and CPCount < 1:
      CMax = npArray[maxInd[i]]
      CPCount+=1
  BMin = -999
  EOSMin = -999 
  for i in range(len(minInd)):
    Point = minInd[i]*PixelCon
    if Point >= .03 and Point < .04:
      BMin = npArray[minInd[i]]
    elif Point >= .075 and Point <= .09:
      EOSMin = npArray[minInd[i]]

  Asy = BMin-EOSMin
  Output = [BMin,EOSMin,CMax,Asy]
  return Output
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#Main Routine
#------------------------------------------------------------------------------

def main():
  """
  The main routine
  """
  ROOT.gStyle.SetOptStat(0)
  if sys.version_info[0] >= 3:
    print ("ERROR: PyROOT only works for Python 2.x")
    raise Exception(" Python version too high")

  #Load in file path name
  nargv = len(sys.argv)
  if nargv < 3:
    inputCSV = "extras/SliceFinder_Example_Files/FEATemps.csv"
  else:
    inputCSV = sys.argv[2] 

  inputfile = sys.argv[1]

  #Make Filename
  try:
    if 'Rec-' in inputfile:
      strFileNum = inputfile.split('Rec-')[1]
      strFileNum = strFileNum[:6] #This gets the number after the recording
      name = 'Rec'+strFileNum
    elif 'plot-' in inputfile:
      strFileNum = inputfile.split('plot-')[1]
      name = 'Rec'+strFileNum
    elif 'Example' in inputfile:
      name = 'Example'
  except:
    name = 'OtherData'

  #Make outdir in file location
  fileloc = inputfile.split('/')[-1]
  fileloc = inputfile.replace(fileloc,'')
  outdir = fileloc + 'ModulesProf/'
  if not os.path.isdir(outdir):
    os.mkdir(outdir)

  c0 = ROOT.TCanvas("c0")
  c0.cd()

  #Get Stave and Data
  Stave = LoadStave(inputfile,name,outdir)
  Slices = GetSlices(Stave,13,outdir)
  GivnSlices = LoadSlices(inputCSV,outdir)

  RInfo = [0 for i in range(13)]
  FInfo = [0 for i in range(14)]


  c0 = ROOT.TCanvas("individual plots","",1400,1000)
  c1 = ROOT.TCanvas("combined plots","",2100,500)
  c1.Divide(3,1)
  TList = ROOT.TList() 
 
  #Create Comparison Plots
  for i in range(GivnSlices[1]):
    RData = Slices[0]
    Rny = Slices[2] 
    RPixCon = Slices[4]
    FData = GivnSlices[0]
    Fny = GivnSlices[2]
    FPixCon = GivnSlices[4]

    #Prints the last real slice twice
    j = i
    if i == GivnSlices[1]-1:
      j = i-1
    else:
      RInfo[j] = FindMinMaxAsy(RData[j][:],RPixCon)
    FInfo[i] = FindMinMaxAsy(FData[i][:],FPixCon)
    
    histR = ROOT.TH1F("histR"+str(j+1),"Module "+str(j+1)+";Position [mm]; Temperature [C]",Rny,0,Rny*RPixCon*1000)
    histF = ROOT.TH1F("histF"+str(i+1),"Simulated Module "+str(i+1)+";Position [mm]; Temperature [C]",Fny,0,Fny*FPixCon*1000)
    #To save them to memory
    TList.Add(histR)
    TList.Add(histF)
    
    histR.SetLineColor(1)
    histR.SetLineWidth(3)
    histF.SetLineColor(2)
    histF.SetLineWidth(3)
   
    #Fill the Histograms
    for y in range(Rny):
      histR.SetBinContent(y+1,RData[j][y])
    for y in range(Fny):
      histF.SetBinContent(y+1,FData[i][y])
    
    #Print individual plots
    c0.cd()
    histR.Draw()
    histR.GetYaxis().SetRangeUser(-40,0)
    if j > 0:
      histR.GetXaxis().SetRangeUser(0,116)

    histF.Draw("Same")
    legend = ROOT.TLegend(0.1,0.7,0.45,0.9)
    legend.AddEntry(histR,"Meas Module "+str(j+1),"l")
    legend.AddEntry(histF,"Sim Module "+str(i+1),"l")
    legend.Draw() 
    c0.Print(outdir+"Comparison"+str(i+1)+".png") 

    #Print triplet graph
    frame = i % 3
    c1.cd(frame+1)
    histR.Draw()
    histR.GetYaxis().SetRangeUser(-40,0)
    if j > 0:
      histR.GetXaxis().SetRangeUser(0,116)

    histF.Draw("Same")
    legend = ROOT.TLegend(0.1,0.7,0.45,0.9)
    legend.AddEntry(histR,"Meas","l")
    legend.AddEntry(histF,"Sim","l")
    legend.Draw() 

    if frame == 2:
      c1.Update()
      c1.Print(outdir + "Comparisons"+str(i-1)+'-'+str(i+1)+".png")

  #Print changes over the stave
  RMods = np.array([i for i in range(1,14)],dtype=float)
  FMods = np.array([i for i in range(1,15)],dtype=float)
 
  RBMin = np.array([RInfo[i][0] for i in range(len(RInfo))])
  FBMin = np.array([FInfo[i][0] for i in range(len(FInfo))])
  REOSMin = np.array([RInfo[i][1] for i in range(len(RInfo))])
  FEOSMin = np.array([FInfo[i][1] for i in range(len(FInfo))])
  RMax = np.array([RInfo[i][2] for i in range(len(RInfo))])
  FMax = np.array([FInfo[i][2] for i in range(len(FInfo))])
  RAsy = np.array([RInfo[i][3] for i in range(len(RInfo))])
  FAsy = np.array([FInfo[i][3] for i in range(len(FInfo))])

  c2 = ROOT.TCanvas("c2","",1400,1000)
  RBMinG = ROOT.TGraph(13,RMods,RBMin)
  FBMinG = ROOT.TGraph(14,FMods,FBMin)
  REOSMinG = ROOT.TGraph(13,RMods,REOSMin)
  FEOSMinG = ROOT.TGraph(14,FMods,FEOSMin)

  RBMinG.SetMarkerStyle(20)
  RBMinG.SetMarkerColor(1)
  RBMinG.SetTitle("Meas Bot Pipe")
  REOSMinG.SetMarkerStyle(22)
  REOSMinG.SetMarkerColor(1)
  REOSMinG.SetTitle("Meas EOS Pipe")
  FBMinG.SetMarkerStyle(20)
  FBMinG.SetMarkerColor(2)
  FBMinG.SetTitle("Sim Bot Pipe")
  FEOSMinG.SetMarkerStyle(22)
  FEOSMinG.SetMarkerColor(2)
  FEOSMinG.SetTitle("Sim EOS Pipe")

  MinGraph = ROOT.TMultiGraph()
  MinGraph.Add(RBMinG,"APL")
  MinGraph.Add(FBMinG,"APL")
  MinGraph.Add(REOSMinG,"APL")
  MinGraph.Add(FEOSMinG,"APL")

  c2.cd()
  c2.SetGrid(1)
  MinGraph.Draw("AP")
  MinGraph.SetTitle("Pipe Minimums")
  MinGraph.GetXaxis().SetTitle("Module Number")
  MinGraph.GetYaxis().SetTitle("Temperature [#circC]")
  MinGraph.GetYaxis().SetRangeUser(-40,-30)
  c2.BuildLegend()
  c2.Update()
  c2.Print(outdir + "MinimaGraph.png")

  RMaxG = ROOT.TGraph(13,RMods,RMax)
  FMaxG = ROOT.TGraph(14,FMods,FMax)

  RMaxG.SetMarkerStyle(20)
  RMaxG.SetMarkerColor(1)
  RMaxG.SetTitle("Meas Cent Max")
  FMaxG.SetMarkerStyle(20)
  FMaxG.SetMarkerColor(2)
  FMaxG.SetTitle("Sim Cent Max")

  MaxGraph = ROOT.TMultiGraph()
  MaxGraph.Add(RMaxG,"APL")
  MaxGraph.Add(FMaxG,"APL")
  MaxGraph.Draw("AP")
  MaxGraph.SetTitle("Stave Central Maximums")
  MaxGraph.GetXaxis().SetTitle("Module Number")
  MaxGraph.GetYaxis().SetTitle("Temperature [#circC]")
  MaxGraph.GetYaxis().SetRangeUser(-40,-25)
  c2.BuildLegend()
  c2.Update
  c2.Print(outdir + "MaximaGraph.png")

  RAsyG = ROOT.TGraph(13,RMods,RAsy)
  FAsyG = ROOT.TGraph(14,FMods,FAsy)

  RAsyG.SetMarkerStyle(20)
  RAsyG.SetMarkerColor(1)
  RAsyG.SetTitle("Meas Asymmetry")
  FAsyG.SetMarkerStyle(20)
  FAsyG.SetMarkerColor(2)
  FAsyG.SetTitle("Sim Asymmetry")

  AsyGraph = ROOT.TMultiGraph()
  AsyGraph.Add(RAsyG,"APL")
  AsyGraph.Add(FAsyG,"APL")
  AsyGraph.Draw("AP")
  AsyGraph.SetTitle("Stave Pipe Asymmetry ")
  AsyGraph.GetXaxis().SetTitle("Module Number")
  AsyGraph.GetYaxis().SetTitle("Temperature [C#circ]")
  AsyGraph.GetYaxis().SetRangeUser(-4,4)
  c2.BuildLegend()
  c2.Update
  c2.Print(outdir + "AsymmetryGraph.png")

if __name__ == '__main__':
  main()
