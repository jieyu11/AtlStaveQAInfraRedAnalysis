#!/usr/bin/env python3

"""
@run:
  ./defectFinder.py [option] [files]

  [option]:   This can be -v or --verbose, prints all fit results
  [files]:    This can be any number of result.root files. If exactly two
              are specified it will do a comparison of the two images. Otherwise
              it will only plot the first one, and then you can then 
              use user commands to plot the rest.
              
              The output files are saved to either flawplots/ or flawplots/fits/
              which it will create. It gives each plot a name using the name of the
              folder containing the result.root file.

@brief:
  The code reads the output root file produced from frameanal.py to find flaws 
  in the cooling pipe temperature profile. It will only find flaws from around
  1 cm in length to 8cm. Flaws smaller than 1 cm are too small to be reliably
  detected from the background and flaws larger than 8cm give global effects
  that are not taken into account by this flaw finding algorithm (And they are
  generally quite obvious in the spectrum).
  

  Once all flaws are found it produces a plot with all of the defects drawn on
  it, and it prints a table that includes each flaw's, center, width and height.

@notes:
  This code assumes using python 2.7.10
  It will not work for python 3.x, due to lack of ROOT support

@email: wheidorn@iastate.edu
"""

import sys          # system 
import os           # operating system
import ROOT         # ROOT from CERN
import numpy as np  # number in python

from defectFinderToolBox import *
from defectFinderPeakFinder import *

#------------------------------------------------------------------------------
#Text Commands-----------------------------------------------------------------
def ComInfo():
  """
  Prints all of the commands
  """
  print("\n-----Stave Defect Finder-----")
  print("h = print list of commands")
  print("q = quit program")
  print("f = shows all attached files")
  print("Defs = prints all defects from every file")
  print("SpDef(i) = prints ith inputfile")
  print("HnC(i,j) = prints hot and cold comparison between ith and jth inputfiles")
  print("TDiff(i,j,type,poff) = prints difference between ith and jth inputfiles type")
  print("                       where the type is given by temp, mean, width, chi2")
  print("                       and poff is an offset given to the first dataset")       
  print("TDiffScale(i,j,type,poff) = same as above but scales file j to i")
  print("SpRMS(i) = prints out a plot of the spectra fit with a line then subtracted")
  print("OneLine(i,j) = prints temperature along whole ith input file with a sep")
  print("                       j in cm")
  print("OLDiff(i,j,l,s) = prints difference between ith and jth inputfiles")
  print("                       along the whole pipelength. The separation is")
  print("                       given by l in cm, and if S is used it rescales j") 
  print("OLDef(i) = prints the temperature and all flaws on one line")
  print("OLMulti(i,j,...) = prints all of the temp profiles together")




def TextCommands(inputfile,outdirs,fitdirs,canvas,bolVerb):
  """
  This allows the user to give commands to replot items and keep the program and
  to edit what is plotted
  """
  ComInfo() 
  bolQuit = False
  while bolQuit == False:
    Vin = str(input("\nType the command you wish to do:\n"))
    #canvas.Clear()
    #canvas.Update()
    if Vin == 'q':
      bolQuit = True
    elif Vin == 'f':
      for i in range(len(inputfile)):
        print(str(i)+" "+inputfile[i])
    elif Vin == 'Defs':
      #Cnew = ROOT.TCanvas("cnew","cnew",0,0,2000,600)
      Defs = []
      counter = 0
      for i in inputfile:
        Defs = np.append(Defs,GetDefects(i,outdirs[counter],fitdirs[counter],canvas,bolVerb))
        if bolVerb == True:
          DefectAnalysis(Defs,outdirs[0],canvas)
        counter+=1
      #DefectAnalysis(Defs,outdirs[0],canvas)
    elif Vin == 'h':
      ComInfo()
    elif 'SpDef(' in Vin:
      Vin = Vin.lstrip('SpDef(')
      Vin = Vin.rstrip(')')
      try:
        Vin = int(Vin)
        #Cnew1 = ROOT.TCanvas("cnew1","cnew1",0,0,2000,600)
        Def = GetDefects(inputfile[Vin],outdirs[Vin],fitdirs[Vin],canvas,bolVerb)
      except:
        print("NOT Correct format") 
    elif 'HnC(' in Vin:
      Vin = Vin.lstrip('HnC(')
      Vin = Vin.rstrip(')')
      Vin = Vin.split(',')
      #Cnew2 = ROOT.TCanvas("cnew2","cnew2",0,0,2000,1200)
      F = HnCComp(inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdirs[int(Vin[0])],fitdirs[int(Vin[0])],canvas,bolVerb) 

    elif 'TDiff(' in Vin:
      Vin = Vin.lstrip('TDiff(')
      Vin = Vin.rstrip(')')
      Vin = Vin.split(',')
      ninputs = len(Vin)

      Ctop = ROOT.TCanvas("ctop","ctop",0,0,2000,490)
      Cbot = ROOT.TCanvas("cbot","cbot",0,510,2000,490)
      if ninputs == 2:
        try:
          F = PlotDiff(0,inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdirs[int(Vin[0])],Ctop,Cbot,'temperature;1')
        except:
          print("2 Inputs, NOT Correct format")
      elif ninputs == 3:
        try:
          options = ["temperature;1","mean;1","width;1","chi2;1"]
          if 'temp' in Vin[2]: choice = options[0]
          elif 'mean' in Vin[2]: choice = options[1]
          elif 'width' in Vin[2]: choice = options[2]
          elif 'chi2' in Vin[2]: choice = options[3]
          else: choice = options[0]
          F = PlotDiff(0,inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdirs[int(Vin[0])],Ctop,Cbot,choice)
        except:
          print("3 Inputs NOT Correct format")
      elif ninputs == 4:
        try:        
          options = ["temperature;1","mean;1","width;1","chi2;1"]
          if 'temp' in Vin[2]: choice = options[0]
          elif 'mean' in Vin[2]: choice = options[1]
          elif 'width' in Vin[2]: choice = options[2]
          elif 'chi2' in Vin[2]: choice = options[3]
          else: choice = options[0]
          F = PlotDiff(int(Vin[3]),inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdirs[int(Vin[0])],Ctop,Cbot,choice)
        except:
          print("4 Inputs NOT Correct format")

    elif 'TDiffScale(' in Vin:
      Vin = Vin.lstrip('TDiffScale(')
      Vin = Vin.rstrip(')')
      Vin = Vin.split(',')
      ninputs = len(Vin)

      Ctop = ROOT.TCanvas("ctop","ctop",0,0,2000,490)
      Cbot = ROOT.TCanvas("cbot","cbot",0,510,2000,490)
      if ninputs == 2:
        try:
          F = PlotDiffScale(0,inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdirs[int(Vin[0])],Ctop,Cbot,'temperature;1')
        except:
         print("2 Inputs, NOT Correct format")
      elif ninputs == 3:
        try:
          options = ["temperature;1","mean;1","width;1","chi2;1"]
          if 'temp' in Vin[2]: choice = options[0]
          elif 'mean' in Vin[2]: choice = options[1]
          elif 'width' in Vin[2]: choice = options[2]
          elif 'chi2' in Vin[2]: choice = options[3]
          else: choice = options[0]
          F = PlotDiffScale(0,inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdirs[int(Vin[0])],Ctop,Cbot,choice)
        except:
          print("3 Inputs NOT Correct format")
      elif ninputs == 4:
        try:        
          options = ["temperature;1","mean;1","width;1","chi2;1"]
          if 'temp' in Vin[2]: choice = options[0]
          elif 'mean' in Vin[2]: choice = options[1]
          elif 'width' in Vin[2]: choice = options[2]
          elif 'chi2' in Vin[2]: choice = options[3]
          else: choice = options[0]
          F = PlotDiffScale(int(Vin[3]),inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdirs[int(Vin[0])],Ctop,Cbot,choice)
        except:
          print("4 Inputs NOT Correct format")
    elif 'TDiffs(' in Vin:
      Vin = Vin.lstrip('TDiffs(')
      Vin = Vin.rstrip(')')
      Vin = Vin.split(',')
      ninputs = len(Vin)
      options = ["temperature;1","mean;1","width;1","chi2;1"]
      Ctop = ROOT.TCanvas("ctop","ctop",0,0,2000,490)
      Cbot = ROOT.TCanvas("cbot","cbot",0,510,2000,490)
      if ninputs == 3:
        for i in range(len(options)):
            
          F = PlotDiff(int(Vin[2]),inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdirs[int(Vin[0])],Ctop,Cbot,options[i])          
      else:
        print("Not correct number of parameters")
    elif 'SpRMS(' in Vin:
      Vin = Vin.lstrip('SpRMS(')
      Vin = Vin.rstrip(')')
      try:
        Vin = int(Vin)
      except:
        continue
      #CRMS = ROOT.TCanvas("crms","crms",0,0,2000,500)
      F = GetBothRMS(inputfile[Vin],outdirs[Vin],canvas)
       
    elif 'OneLine(' in Vin:
      Vin = Vin.lstrip('OneLine(')
      Vin = Vin.rstrip(')')
      VinList = Vin.split(',')
      try:
        Vin = int(VinList[0])
        Length = float(VinList[1])
      except:
        continue
      COL = ROOT.TCanvas("canvas 1 line","canvas 1 line",0,0,2000,600)
      OneLine(inputfile[Vin],outdirs[Vin],COL,Length)
    elif 'OLDiff(' in Vin:
      Vin = Vin.lstrip('OneLineDiff(')
      Vin = Vin.rstrip(')')
      VinList = Vin.split(',')
      try:
        File1 = int(VinList[0])
        File2 = int(VinList[1])
        try:
          Length = float(VinList[2])
        except:
          Length = 10.
      except:
        continue
      Scale = False
      if len(VinList) > 3 and VinList[3] == 'S':
        Scale = True
      #CLDiff = ROOT.TCanvas("cline","cline",0,0,2000,700)
      OneLineComp(inputfile[File1],inputfile[File2],outdirs[File1],canvas,Length,Scale)
    elif 'OLDef(' in Vin:
      Vin = Vin.lstrip('OneLineDef(')
      Vin = Vin.rstrip(')')
      try:
        File = int(Vin) 
      except:
        continue
      Defs = GetOneLineDefects(inputfile[File],outdirs[File],fitdirs[File],canvas,bolVerb)
    elif 'OLRMS(' in Vin:
      Vin = Vin.lstrip('OneLineRMS(')
      Vin = Vin.rstrip(')')
      try:
        File = int(Vin)
      except:
        continue
      #COLRMS = ROOT.TCanvas("clineRMS","clineRMS",0,0,800,800)
      OneLineRMS(inputfile[File],outdirs[File],canvas)
    elif 'OLMulti(' in Vin:
      Vin = Vin.lstrip('OLMulti(')
      Vin = Vin.rstrip(')')
      VinList = Vin.split(',')
      try:
        Files = []
        Outdirs=[]
        if Vin == 'a':
          Files = inputfile

        else:
          for i in VinList:
            Files = np.append(Files,inputfile[int(i)])
            Outdirs = np.append(Outdirs,outdirs[int(i)])
      except:
        print("Wrong input variables")
        continue
      OneLineMulti(Files,Outdirs,canvas)

    elif 'TempMean(' in Vin: 
      Vin = Vin.lstrip('TempMean(')
      Vin = str(Vin.rstrip(')'))
      VinList = Vin.split(',')
      try:
        Files = [] 
        Outdirs = []
        if Vin == 'A':
          Files = inputfile
        else:
          for i in VinList:
            Files = np.append(Files,inputfile[int(i)])
            Outdirs = np.append(Outdirs,outdirs[int(i)])
      except:
        print("Wrong input variables")
        continue 
      #CTM = ROOT.TCanvas("CTM","cTempMean",0,0,2000,2000)
      counter = 0
      for fyle in Files:
        FindTempHist(str(fyle),Outdirs[counter],canvas)
        counter+=1

    elif 'FindAvgTemps' in Vin:
      #canvas = ROOT.TCanvas("genericCanvas")
      FindAvgTemps(inputfile,outdir,canvas) 
    else:
      print("Command not recognized")

#------------------------------------------------------------------------------
#Main Routine------------------------------------------------------------------
def main():
  """
  The main routine
  """
  
  #Load the Files
  nargv = len(sys.argv)
  inputfile = []
  if (nargv <= 1):
    print("ERROR: Please provide an input root file or set of files. These files should be result.root from frameanal.py")
    return
  elif (nargv >= 3) and (sys.argv[1] == '-v' or sys.argv[1] == '--verbose'):
    print("VERBOSE: Creating a plot for each defect fit")
    bolVerb = 1    
    for i in range(2,nargv):
      inputfile = np.append(inputfile,sys.argv[i])
  else:
    for i in range(1,nargv):
      inputfile = np.append(inputfile,sys.argv[i])
    bolVerb = 0

  #MAKING THE OUTDIRECTORYs
  outdirs=[]
  fitdirs=[]
  for fyle in inputfile:
    basicName = fyle.split("/")[-1]
    pathToFile = fyle.replace(basicName,"")
    outdir = pathToFile + "flawplots/"
    fitdir = outdir+"fits/"
    if not os.path.isdir( outdir ):
      os.mkdir( outdir )
    if not os.path.isdir( fitdir ):
      os.mkdir( fitdir )
    outdirs.append(outdir)
    fitdirs.append(fitdir)

  canvas = ROOT.TCanvas("main canvas","canvas",0,0,2000,600)
  canvas.SetGridx()

  #First Run
  if len(inputfile) == 1:
    Def = GetDefects(inputfile[0],outdirs[0],fitdirs[0],canvas,bolVerb)
    if bolVerb > 0:
      DefectAnalysis(Def,outdirs[0],canvas)
  if len(inputfile) == 2:
    Def = HnCComp(inputfile[0],inputfile[1],outdirs[0],fitdirs[0],canvas,bolVerb)
  else:
    OneLineMulti(inputfile,outdirs,canvas)
   
  #Wait for commands
  TextCommands(inputfile,outdirs,fitdirs,canvas,bolVerb)

if __name__ == "__main__":
  main()
