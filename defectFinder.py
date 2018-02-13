#!/usr/bin/python

"""
@run:
  ./defectFinder.py plot/result.py

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


#------------------------------------------------------------------------------
#Some Functions----------------------------------------------------------------
def GetHistogram(inputfile,intPipeNum,strType):
  """
  Gets the histogram from the input file and returns it
  """
  File = ROOT.TFile(inputfile,"READ")
  side = "top_pipe_"
  if intPipeNum == 1:
    side = "bot_pipe_"

  Hist = File.Get(side+strType)
  Hist.SetDirectory(0)
  File.Close()
  return Hist

#------------------------------------------------------------------------------
def GetHistInfo(Histogram):
  """
  Finds the number of bins, start and end points of a histogram
  """
  nbins = Histogram.GetNbinsX()
  Xmin = Histogram.GetXaxis().GetBinUpEdge(0)
  Xmax = Histogram.GetXaxis().GetBinUpEdge(nbins)
  Name = Histogram.GetName()
  return [nbins,Xmin,Xmax,Name]

#------------------------------------------------------------------------------
#FindPeaks---------------------------------------------------------------------

def FindPeaks(files,outdir,fitdir,canvas,intPipeNum,bolVerb = 0,Boixes = [], strType = "temperature;1"):
  """
  Uses the information along the cooling pipe to find the flaws along the pipe.
  It has three main parts.

  1. Histogram conversion:
    - Loads Histogram from file.
    - If necessary, flips it so that there are peaks
    - Band passes the spectrum to remove high and low fluctuations
    - Uses TSpectrum.Search to find peaks
    - Uses TSpectrum.Background to find background histogram
    - Subtracts Background from inverted bandpassed histogram to get HistPeaks
  2. Fitting to find flaws:
    At each peak value found by TSpectrum the following occurs
    - Fit 1: A fit to HistPeaks using a gaussian function
    - Fit 2: A fit to HistPeaks using the width and height found from Fit 1 
              using a gaussian with an offset.
    - Fit 3: A fit to HistInverted using the width and height from the best of
              Fit 1 and 2. This fit is only kept if it is better than the best
              previous fit.
    - Height Finding: Fit 3 is used again on HistInverted2
             (HistInverted without the bandpass filter). 
    - Cuts:
        @ Position- TSpectrum peak value must be within 2 cm of fitted value
                     Peak maximum must not be within 2 cm of edge
        @ Peak Width must be > 1cm and < 8cm
        @ Peak Height must be > 0.2 C
  3. Plotting and output:
    Each peak that passes the cut is then saved as a defect and plotted as a
    box on the input spectrum. The flaw information is then returned by the
    function.
  """
  strPipeLab = "top_"
  if intPipeNum == 1:
    strPipeLab = "bot_"
  Hist=GetHistogram(files,intPipeNum,strType)
  Info = GetHistInfo(Hist)
  nBins = Info[0] 
  Xmin = Info[1] 
  Xmax =  Info[2] 
  npeaks = 20
  objSpectrum = ROOT.TSpectrum(2*npeaks)
  HistInverted =  Hist.Clone()
     
  for i in range(nBins):
    HistInverted.SetBinContent(i+1,-1*Hist.GetBinContent(i+1)) 
  HistInverted2 = HistInverted.Clone() 

  #FINDING FLAWS---------------------------------------------------------- 
  HistInverted.Draw()
  HistInverted = BandPassFFT(HistInverted,2,200) 
  HistBackground = objSpectrum.Background(HistInverted,17,"")

  HistPeaks = HistInverted - HistBackground

  sigma = 6
  threshold = 0.05
  nfound = objSpectrum.Search(HistInverted,sigma,"noMarkov",threshold)
  peakPosX = objSpectrum.GetPositionX()

  HistInverted.Draw()
  HistBackground.Draw("same")
  if bolVerb > 0:
    canvas.Print(fitdir+files[5:11]+strPipeLab+"BackgroundInvertPlot.png")
 
  PeakPosArray = np.zeros(nfound) 
  for i in range(nfound):
    PeakPosArray[i] = peakPosX[i]  
  #print (str(PeakPosArray))
  
  MFlaws = []
  for i in range(nfound):
    #Find the Flaws Method 1---------- 
    fit = ROOT.TF1("pgaus1","[0]*exp(-0.5*((x-[1])/[2])^2)",Xmin,Xmax)
    fit.SetParLimits(0,0,2)
    fit.SetParLimits(1,peakPosX[i]-5,peakPosX[i]+5)
    fit.SetParLimits(2,0.25,10)
    fit.SetParameters(1,peakPosX[i],1) 
    HistPeaks.Fit(fit,"Q","",peakPosX[i]-1.5*sigma,peakPosX[i]+1.5*sigma)
    #Get the fit info
    FitPeakPosX = fit.GetParameter(1)
    FitSigma = abs(fit.GetParameter(2))
    FitChiSq = fit.GetChisquare()
    FitDegFr = fit.GetNDF()
    FitHeight = abs(fit.GetParameter(0)) 
    FitGoodness = abs(FitChiSq/FitDegFr)
    FitLevel = 0

    #Find the Flaws Method 2----------
    fit2 = ROOT.TF1("pgaus2","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",Xmin,Xmax)
    fit2.SetParLimits(0,0,5)
    fit2.SetParLimits(1,peakPosX[i]-5,peakPosX[i]+5)
    fit2.SetParLimits(2,0.25,10)
    fit2.SetParLimits(3,-.5,.5)
    fit2.SetParameters(FitHeight,FitPeakPosX,FitSigma,0) 
    HistPeaks.Fit(fit2,"Q","",peakPosX[i]-sigma,peakPosX[i]+sigma)
    #Get the fit info
    Fit2PeakPosX = fit2.GetParameter(1)
    Fit2Sigma = abs(fit2.GetParameter(2))
    Fit2ChiSq = fit2.GetChisquare()
    Fit2DegFr = fit2.GetNDF()
    Fit2Height = abs(fit2.GetParameter(0)) 
    Fit2Goodness = abs(Fit2ChiSq/Fit2DegFr)

    if Fit2Goodness < FitGoodness:
      FitPeakPosX = Fit2PeakPosX
      FitSigma = Fit2Sigma
      FitHeight = Fit2Height
      FitGoodness = Fit2Goodness
      FitLevel = 1
      if bolVerb > 0:
        canvas.Print(fitdir+files[5:11]+strPipeLab+str(i)+"_FitResults.png")

    #Find Flaws Method 3----------
    fit3 = ROOT.TF1("pgaus3","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",Xmin,Xmax)
    fit3.SetParLimits(0,FitHeight-0.1,FitHeight+0.1)
    fit3.SetParLimits(1,FitPeakPosX-1,FitPeakPosX+1)
    fit3.SetParLimits(2,FitSigma-1,FitSigma+1)
    fit3.SetParLimits(3,-1,1)
    fit3.SetParLimits(4,-2,2)
    fit3.SetParameters(FitHeight,FitPeakPosX,FitSigma,0,0) 
    HistInverted.Fit(fit3,"Q","",peakPosX[i]-sigma,peakPosX[i]+sigma)
    #Get the fit info
    Fit3PeakPosX = fit3.GetParameter(1)
    Fit3Sigma = abs(fit3.GetParameter(2))
    Fit3ChiSq = fit3.GetChisquare()
    Fit3DegFr = fit3.GetNDF()
    Fit3Height = abs(fit3.GetParameter(0)) 
    Fit3Goodness = abs(Fit3ChiSq/Fit2DegFr)

    if Fit3Goodness < FitGoodness:
      FitPeakPosX = Fit3PeakPosX
      FitSigma = Fit3Sigma
      FitHeight = Fit3Height
      FitGoodness = Fit3Goodness
      FitLevel = 2
      if bolVerb > 0:
        canvas.Print(fitdir+files[5:11]+strPipeLab+str(i)+"_FitResults.png")

    #Height Finding
    fit4 = ROOT.TF1("pgaus4","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",Xmin,Xmax)
    fit4.SetParLimits(0,FitHeight-1,FitHeight+2)
    fit4.SetParLimits(1,FitPeakPosX-1,FitPeakPosX+1)
    fit4.SetParLimits(2,FitSigma-0.5,FitSigma+0.5)
    fit4.SetParLimits(3,-50,50)
    fit4.SetParLimits(4,-2,2)
    fit4.SetParameters(FitHeight,FitPeakPosX,FitSigma,0,0) 
    HistInverted2.Fit(fit4,"Q","",peakPosX[i]-sigma,peakPosX[i]+sigma)
    FitHeightMeas = fit4.GetParameter(0)

    AcceptX = 2         #   Any peak that is found that is not +/- acceptace from
                        # the edges of the histogram or the peakX position 
                        # found by TSpectrum are discarded
    SizeConstant = 2.5  #   To be at 95% of the peaks area, one uses 4 
                        # (This seems to be about right for the flaw size 
                        # to implimented size ratio)

    if FitPeakPosX > peakPosX[i] -AcceptX and FitPeakPosX < peakPosX[i]+AcceptX and\
       FitPeakPosX>Xmin+AcceptX and FitPeakPosX< Xmax-AcceptX and\
       FitGoodness < 0.05:
      if abs(FitSigma*SizeConstant)<8 and abs(FitSigma*SizeConstant) > 1 and FitHeight > 0.2: 
        MFlaws = np.append(MFlaws,[FitLevel,FitPeakPosX,FitSigma,FitHeightMeas,FitGoodness])
  
  #PLOTTING---------------------------------------------------------------
  fitline = ROOT.TF1("offset","[0]",Xmin,Xmax)
  Hist.Fit(fitline,"Q N","",Xmin,Xmax)
  Offset = fitline.GetParameter(0)

  canvas.Update()
  Hist.Draw()

  #Histogram Y Range
  if Offset > 22:
    Hist.SetMaximum(Offset +1)
    Hist.SetMinimum(Offset -3)
  elif Offset < -10:
    Hist.SetMaximum(Offset +4)
    Hist.SetMinimum(Offset -1)
  elif Offset < 10:
    Hist.SetMaximum(Offset +3)
    Hist.SetMinimum(Offset -1)
  else:
    Hist.SetMaximum(Offset+2)
    Hist.SetMinimum(Offset-2)
    

  nMflaws = len(MFlaws)/5
 
  global Boxes 
  BoxesL = len(Boxes)
 
  Ymax = Hist.GetMaximum()
  Ymin = Hist.GetMinimum()

  #Major Flaw Boxes
  for i in range(nMflaws):
    Center = MFlaws[i*5+1]
    Xmin = Center - SizeConstant/2*MFlaws[i*5+2]
    Xmax = Center + SizeConstant/2*MFlaws[i*5+2] 
    Boxes = np.append(Boxes,ROOT.TBox(Xmin,Ymin,Xmax,Ymax))
    Boxes[i+BoxesL].SetLineColor(i+2)
    Boxes[i+BoxesL].SetLineWidth(1)
    FitGoodness = MFlaws[i*5+4]
    trans = 0.5
    if FitGoodness < 0.005:
      Boxes[i+BoxesL].SetFillColorAlpha(2,trans)
    elif FitGoodness < 0.01:
      Boxes[i+BoxesL].SetFillColorAlpha(5,trans)
    elif FitGoodness < 0.03:
      Boxes[i+BoxesL].SetFillColorAlpha(3,trans)
    else:
      Boxes[i+BoxesL].SetFillColorAlpha(4,trans)
    Boxes[i+BoxesL].Draw()

  if bolVerb > 0: 
    canvas.Print(outdir+files[5:11]+strPipeLab+"Tpeaks.root")
    canvas.Print(outdir+files[5:11]+strPipeLab+"Tpeaks.png")

  #PREPPINGDEFECTINFO-----------------------------------------------------
  DefectInfo = []
  for i in range(nMflaws):
    DefectInfo = np.append(DefectInfo,i)                                  #Flaw Number
    DefectInfo = np.append(DefectInfo,intPipeNum)                         #Flaw Pipe
    DefectInfo = np.append(DefectInfo,float(MFlaws[i*5+1]))               #Flaw Center
    DefectInfo = np.append(DefectInfo,float(SizeConstant*MFlaws[i*5+2]))  #Flaw Width
    DefectInfo = np.append(DefectInfo,float(MFlaws[i*5+3]))               #Flaw Height
    DefectInfo = np.append(DefectInfo,float(MFlaws[i*5+4]))               #Fit Goodness
    DefectInfo = np.append(DefectInfo,MFlaws[i*5])                        #Fit Level 
  return DefectInfo

#------------------------------------------------------------------------------
def BandPassFFT(Hist,CutOffL,CutOffH):
  """
    This filter switches the temperature data to frequency space via the hfft 
    method. The high and low values given by the users replace all the frequency
    values outside the range with zeros. It is then inversed back into real space
    and the histogram is given back.
  """
  Info = GetHistInfo(Hist)
  nbins = Info[0]
  Xmin = Info[1]
  Xmax = Info[2]
  Data = np.zeros(nbins)
  for i in range(nbins):
    Data[i] = Hist.GetBinContent(i+1)
  FFTData = np.real(np.fft.hfft(Data))

  if CutOffL < 0:
    CutOffL = 0
  if CutOffH > 2*nbins-2:
    CutOffH = 2*nbins-2

  for i in range(2*nbins-2):
    if i > CutOffH:
      FFTData[i] = 0.0
    if i < CutOffL:
      FFTData[i] = 0.0
  
  FiltData = np.fft.ihfft(FFTData)
  HistHPass = ROOT.TH1F(Info[3]+"FFT",Info[3]+"FFT",nbins,Xmin,Xmax)
  for i in range(nbins):
    HistHPass.AddBinContent(i+1,np.real(FiltData[i]))
  return HistHPass  

#------------------------------------------------------------------------------
def GetDefects(inputfile,outdir,fitdir,C1,bolVerb = 0):
  """
  This takes an input file and gets all of the defects and prints all of the 
  relevant information. This is necessary for printing the boxes on a single canvas
  """
  DefectInfo = []
  C1.Clear()
  C1.Divide(1,2)
  C1.cd(1)
  global Boxes
  Boxes = []
  DefectInfo = np.append(DefectInfo,FindPeaks(inputfile,outdir,fitdir,C1,0,bolVerb,Boxes)) 
  DefectsTop = len(Boxes)
  for i in range(DefectsTop):
    Boxes[i].Draw()
  C1.cd(2)
  DefectInfo = np.append(DefectInfo,FindPeaks(inputfile,outdir,fitdir,C1,1,bolVerb,Boxes))
  for i in range(len(Boxes)-DefectsTop):
    Boxes[i+DefectsTop].Draw()
  C1.Print(outdir+inputfile[5:11]+"AllFlaws.png")
  C1.Print(outdir+inputfile[5:11]+"AllFlaws.root")



  PrintDefectInfo(DefectInfo,-1)
  return DefectInfo

#------------------------------------------------------------------------------
def PrintDefectInfo(DefectInfo,intLine):
  """
  This prints the DefectInfo array out into an easily readable format. 
  If intLine = -1 it prints all the information
  """
  bolPrintAll = False
  if intLine == -1:
    bolPrintAll = True
  print("-------------------Flaws Found-------------------------")
  print(" Name    Center     Width     Height    Goodness    Fit")
  nThings = 7
  if bolPrintAll == True:
    for i in range(len(DefectInfo)/nThings):
      print " {0:2}-{1:1}  {2:8.3f}  {3:8.3f}  {4:8.3f}  {5:8.5f}    {6:1}".format(int(DefectInfo[i*nThings]),int(DefectInfo[i*nThings+1]),DefectInfo[i*nThings+2],DefectInfo[i*nThings+3],DefectInfo[i*nThings+4],DefectInfo[i*nThings+5],int(DefectInfo[i*nThings+6]))
  else:  
    print " {0:2}-{1:1}  {2:8.3f}  {3:8.3f}  {4:8.3f}  {5:8.5f}    {6:1}".format(int(DefectInfo[intLine*nThings]),int(DefectInfo[intLine*nThings+1]),DefectInfo[intLine*nThings+2],DefectInfo[intLine*nThings+3],DefectInfo[intLine*nThings+4],DefectInfo[intLine*nThings+5],int(DefectInfo[intLine*nThings+6]))
  if len(DefectInfo) == 0:
    print(" NO FLAWS WERE FOUND ")
#------------------------------------------------------------------------------
def DefectAnalysis(Defects,outdir,Canvas):
  """
  Takes in a set of defect data, strips out each set to create many histograms and plots
  """
  nThings = 7
  nDefects = len(Defects)/nThings  

  #Make Histograms of Width, Height, Goodness, Fit
  Name = ["Width","Height","Goodness","Fit"]
  Xmin = [0,0,0,0]
  Xmax = [10,5,0.05,3]
  nbins = 20
  Canvas.cd()

  for i in range(4):
    Hist = ROOT.TH1F(Name[i]+"TH1F",Name[i],nbins,Xmin[i],Xmax[i])
    for j in range(nDefects):
      Hist.Fill(Defects[j*nThings+i+3],1)
    Hist.Draw()
    Canvas.Print(outdir+Name[i]+"_Hist.png")

  Width =[]
  Height = []
  Goodness =[]
  for i in range(nDefects):
    Width = np.append(Width,Defects[i*nThings+3])
    Height = np.append(Height,Defects[i*nThings+4])
    Goodness = np.append(Goodness,Defects[i*nThings+5])
  print(Width)
  mg = ROOT.TMultiGraph()
  Canvas.Clear()
  WHGraph = []
  for i in range(nDefects):
    W = []
    H = []
    W = np.append(W,Width[i])
    H = np.append(H,Height[i])
    WHGraph = np.append(WHGraph,ROOT.TGraph(1,W,H))
    WHGraph[i].SetMarkerSize(1.5)
    WHGraph[i].SetMarkerColor(1)
    WHGraph[i].SetMarkerStyle(21)
    #WHGraph[i].Draw("AP same")
    mg.Add(WHGraph[i])


  mg.Draw("AP")
  mg.GetXaxis().SetTitle("Width [ cm ]")
  mg.GetYaxis().SetTitle("Height [ #circC ]")
  
  Canvas.Print(outdir+"WidthVHeightGraph.png")
  Canvas.Print(outdir+"WidthVHeightGraph.root")

#------------------------------------------------------------------------------
def main():
  """
  The main routine
  """
  if sys.version_info[0] >= 3:
    print ("ERROR: PyROOT only works with Python 2.x")
    raise Exception(" Python Version too high. Use 2.x.")

  #Load the Files
  nargv = len(sys.argv)
  if (nargv <= 1):
    print("ERROR: Please provide an input root file. These files should be result.root from frameanal.py")
    return
  elif (nargv >= 3) and (sys.argv[1] == '-v' or sys.argv[1] == '--verbose'):
    print("VERBOSE: Creating a plot for each defect fit")
    bolVerb = 1    
    inputfile = sys.argv[2]
  else:
    bolVerb = 0
    inputfile = sys.argv[1]

  #MAKING THE OUTDIRECTORY
  outdir = "flawplots/"
  fitdir = outdir+"fits/"
  if not os.path.isdir( outdir ):
    os.mkdir( outdir )
  if not os.path.isdir( fitdir ):
    os.mkdir( fitdir )

  C1 = ROOT.TCanvas("c1","c1",0,0,2000,600)
  
  Def = GetDefects(inputfile,outdir,fitdir,C1,bolVerb)

  if bolVerb > 0:
    DefectAnalysis(Def,outdir,C1)

  """
  #-------------------------------
  #These Lines can be used as a module that looks at many results
  #FINDING DEFECTS
  print("Stave2RJ----687--------")
  Def = GetDefects("plot-000687-2RJ-p50/result.root",outdir,fitdir,C1,bolVerb)
  print("Stave2RL----686--------")
  Def = np.append(Def,  GetDefects("plot-000686-2RL-p50/result.root",outdir,fitdir,C1,bolVerb))
  print("Stave5J----018---------")
  Def = np.append(Def,  GetDefects("plot-000018-5J-p50/result.root",outdir,fitdir,C1,bolVerb))
  print("Stave5L----019---------")
  Def = np.append(Def,  GetDefects("plot-000019-5L-p50/result.root",outdir,fitdir,C1,bolVerb))
  print("Stave2J----701---------")
  Def = np.append(Def,  GetDefects("plot-000701-2J-p50/result.root",outdir,fitdir,C1,bolVerb))
  print("Stave2L----703---------")
  Def = np.append(Def,  GetDefects("plot-000703-2L-p50/result.root",outdir,fitdir,C1,bolVerb)) 
  
  DefectAnalysis(Def,outdir,C1)
  #-------------------------------
 """ 


if __name__ == "__main__":
  main()
