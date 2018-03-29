#!/usr/bin/python

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

#------------------------------------------------------------------------------
#Some Functions----------------------------------------------------------------
def GetHistogram(inputfile,intPipeNum,strType = "temperature;1"):
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
def ShiftHistogram(HistToShift,intPixelsToShift):
  """
  Cuts a side of the histogram by a number of pixels 0=Left 1=Right
  """ 
  Info=GetHistInfo(HistToShift)
  C0 = 0.23272727
  nbins = int(Info[0])
  Xmin = float(Info[1] + C0*intPixelsToShift)
  Xmax = float(Info[2] + C0*intPixelsToShift)
  HistShifted = ROOT.TH1F(Info[3]+"Shifted",Info[3]+"Shifted",nbins,Xmin,Xmax)
  for i in range(nbins):
    HistShifted.AddBinContent(i+1,HistToShift.GetBinContent(i+1))
  return HistShifted


#------------------------------------------------------------------------------
def CutHistogram(HistToCut,intSide,intPixelsToCut):
  """
  Cuts a side of the histogram by a number of pixels 0=Left 1=Right
  """ 
  Info=GetHistInfo(HistToCut)
  C0 = 0.23272727
  nbins = int(Info[0]-intPixelsToCut)

  strSide="LCut"
  Xmin = float(Info[1]+ C0*intPixelsToCut)
  Xmax = Info[2]
  if intSide == 1:
    strSide="RCut"
    Xmin = Info[1]
    Xmax = float(Info[2] - C0*intPixelsToCut)

  HistCut = ROOT.TH1F(Info[3]+strSide,Info[3]+strSide,nbins,Xmin,Xmax)
  if intSide == 0:
    for i in range(nbins):
      HistCut.AddBinContent(i+1,HistToCut.GetBinContent(i+intPixelsToCut+1))
    return HistCut
  elif intSide == 1:
    for i in range(nbins):
      HistCut.AddBinContent(i+1,HistToCut.GetBinContent(i+1))
    return HistCut
  
#------------------------------------------------------------------------------
def CutHistToSameSize(Hist1,Hist2):
  """
  cuts the first histogram to the area of the second and outputs the cut histogram
  """
  C0 = 0.23272727
  Info1 = GetHistInfo(Hist1)
  Info2 = GetHistInfo(Hist2)
  nBinDif = abs(Info1[0] - Info2[0])
  XminDif = Info1[1]-Info2[1]
  XmaxDif = Info1[2]-Info2[2]

  #Cut Left Side
  if XminDif > 0:
    Hist2 = CutHistogram(Hist2,0,int(round(XminDif/C0)))
  elif XminDif < 0:
    Hist1 = CutHistogram(Hist1,0,int(round(abs(XminDif)/C0)))

  #Cut Right Side
  if XmaxDif > 0:
    Hist1 = CutHistogram(Hist1,1,int(round(XmaxDif/C0)))
  elif XmaxDif < 0:
    Hist2 = CutHistogram(Hist2,1,int(round(abs(XmaxDif)/C0)))

  return [Hist1,Hist2]

#------------------------------------------------------------------------------
def TempIsHot(Hist):
  """
    This takes a histogram and finds if it is hot or cold 
  """
  Info = GetHistInfo(Hist)
  nBins = Info[0]
  Xmin = Info[1]
  Xmax = Info[2]
  fittemp = ROOT.TF1("Offset","[0]",Xmin,Xmax)
  Hist.Fit(fittemp,"Q nodraw","",Xmin,Xmax)
  AvgTemp = fittemp.GetParameter(0)
  if AvgTemp > 30:
    return True
  elif AvgTemp < 10:
    return False
  else:
    print("WARNING: Temperature values are too close to ambient for decent results")

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
def MakeFileName(inputfile):
  """
  Takes the input file and makes a name for it
  """
  filename = inputfile.split('/')
  try: 
    filename = filename[-2]
    if '.root' in filename:
      raise Exception(".root in output filename")
  except:
    #If it cannot use a file structure it creates a date for the name
    import datetime
    filename = datetime.datetime.now().strftime("%I-%M-%S")
  return filename

#------------------------------------------------------------------------------
#Plotting Difference Between Data Sets-----------------------------------------
def PlotDiffLine(offset,file1,file2,outdir,objCanvas,intPipeNum,strType):
  """
  creates a canvas, plots both histograms and their difference. If an offset value
  is given, it will shift the first histogram that many pixels to the right.
  """
  Hist0 = GetHistogram(file1,intPipeNum,strType) 
  Hist1 = GetHistogram(file2,intPipeNum,strType)
  Hist0info = GetHistInfo(Hist0)
  Hist1info = GetHistInfo(Hist1)

  if offset != 0:
    Hist0 = ShiftHistogram(Hist0,offset)
  
  #Check to see if histograms are not on the same boundaries, fix them if they are.
  while Hist0info[0:2] != Hist1info[0:2]:
    Hists = CutHistToSameSize(Hist0,Hist1)
    Hist0 = Hists[0]
    Hist1 = Hists[1]
    Hist0info = GetHistInfo(Hist0)
    Hist1info = GetHistInfo(Hist1)

  #Get the Difference Plot
  HistDiff = Hist1-Hist0

  nbins = HistDiff.GetNbinsX()
  fit = ROOT.TF1( "FitLine", "pol0", 0., float(nbins) )

  HistDiff.Fit(fit,"q")
  pars = np.zeros(3,dtype=float)
  fit.GetParameters(pars)
  offset = pars[0]
  Xmin= Hist1.GetXaxis().GetBinUpEdge(0)
  Xmax= Hist1.GetXaxis().GetBinUpEdge(nbins)
  HistOffset = ROOT.TH1F("OffsetHistogram","Offset",int(nbins),Xmin,Xmax)
  
  for i in range(nbins):
    HistOffset.AddBinContent(i,offset)

  Hist0Shifted = Hist0
  #Hist1Shifted = Hist1-HistOffset
  Hist1Shifted = Hist1
  
  yMax = max(Hist0Shifted.GetMaximum(),Hist1Shifted.GetMaximum())
  yMin = min(Hist0Shifted.GetMinimum(),Hist1Shifted.GetMinimum())
 
  ySep = yMax-yMin

  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)

  #Pad Basics
  LS = 0.15
  TitleS = 0.15
  TitleO =0.7

  strPipeLab = "top_"
  if intPipeNum == 1:
    strPipeLab = "bot_"
  Title = strPipeLab+strType[:-2]

  Pad1 = ROOT.TPad(strPipeLab+"pad1","Both Plots",0.0,0.3,1.0,1.0)
  Pad2 = ROOT.TPad(strPipeLab+"pad2","Difference Value",0.0,0.0,1.0,0.3)
  Pad1.Draw()
  Pad1.SetGrid()
  Pad2.Draw()
  Pad2.SetGrid()

  Hist0Title = MakeFileName(file1)
  Hist1Title = MakeFileName(file2)

  #Plot 1 Comparison
  Pad1.cd()
  Hist1Shifted.SetLineColor(1)
  Hist0Shifted.Draw("LC CM")
  Hist0Shifted.SetAxisRange(yMin-0.1*ySep,yMax+0.1*ySep,"Y")
  Hist1Shifted.Draw("SAME LC CM")
  Hist1Shifted.SetAxisRange(yMin-0.1*ySep,yMax+0.1*ySep,"Y") 
  Hist0Shifted.GetYaxis().SetLabelSize(LS*.3)
  Hist0Shifted.GetYaxis().SetTitleSize(TitleS*.3)
  Hist0Shifted.GetYaxis().SetTitleOffset(TitleO)
  Hist0Shifted.GetXaxis().SetLabelSize(LS*.3)
  Hist0Shifted.GetXaxis().SetTitleSize(TitleS*.3)
  Hist0Shifted.GetXaxis().SetTitleOffset(TitleO)

  Legend = ROOT.TLegend(0.1,0.7,0.3,0.9)
  Legend.SetHeader(Title,"C")
  Legend.AddEntry(Hist0Shifted,Hist0Title,"L")
  Legend.AddEntry(Hist1Shifted,Hist1Title,"L")
  Legend.Draw() 

  #Plot 1 Delta
  Pad2.cd()
  HistDiff.Draw()
  HistDiff.GetYaxis().SetTitle(strType[:-2]+" Difference")
  HistDiff.GetYaxis().SetLabelSize(LS*.7)
  HistDiff.GetYaxis().SetTitleSize(LS*.7) 
  HistDiff.GetYaxis().SetTitleOffset(TitleO*0.5)
  #HistDiff.SetAxisRange(offset-2,offset+2,'Y')
  HistDiff.GetXaxis().Delete()
  HistDiff.SetNdivisions(10,"Y")
  HistDiff.GetXaxis().SetLabelOffset(10)
  filename = MakeFileName(file1)+MakeFileName(file2)
  
  #PRINTING THE PLOTS  
  objCanvas.Print(outdir+filename+strPipeLab+strType[:-2]+"_dif.png")  
  objCanvas.Print(outdir+filename+strPipeLab+strType[:-2]+"_dif.root")
 

#------------------------------------------------------------------------------
def PlotDiff(offset,file1,file2,outdir,C1,C2,strType = "temperature;1"):
  """
  Plots the difference between two stave image temperatures along both the top and
  bottom cooling pipe.
  """

  C1.cd()
  PlotDiffLine(offset,file1,file2,outdir,C1,0,strType)
  
  C2.cd() 
  PlotDiffLine(offset,file1,file2,outdir,C2,1,strType) 

#------------------------------------------------------------------------------
#FindPeaks---------------------------------------------------------------------

def FindPeaks(files,outdir,fitdir,canvas,intPipeNum,bolVerb = 0, strType = "temperature;1"):
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
    - Fit 3: A fit to HistClone using the width and height from the best of
              Fit 1 and 2. This fit is only kept if it is better than the best
              previous fit.
    - Height Finding: Fit 3 is used again on HistClone2
             (HistClone without the bandpass filter). 
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
  HistClone =  Hist.Clone()

  filename = MakeFileName(files)
 
  bolTempIsHot = TempIsHot(Hist)
  if bolTempIsHot == True:
    for i in range(nBins):
      HistClone.SetBinContent(i+1,-1*Hist.GetBinContent(i+1)) 
  HistClone2 = HistClone.Clone() 

  #FINDING FLAWS---------------------------------------------------------- 
  HistClone.Draw()
  HistClone = BandPassFFT(HistClone,2,200) 
  HistBackground = objSpectrum.Background(HistClone,17,"")

  HistPeaks = HistClone - HistBackground

  sigma = 6
  threshold = 0.05
  nfound = objSpectrum.Search(HistClone,sigma,"noMarkov",threshold)
  peakPosX = objSpectrum.GetPositionX()

  HistClone.Draw()
  HistBackground.Draw("same")
  if bolVerb > 0:
    canvas.Print(fitdir+filename+strPipeLab+"BackgroundInvertPlot.png")
 
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
        canvas.Print(fitdir+filename+strPipeLab+str(i)+"_FitResults.png")

    #Find Flaws Method 3----------
    fit3 = ROOT.TF1("pgaus3","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",Xmin,Xmax)
    fit3.SetParLimits(0,FitHeight-0.1,FitHeight+0.1)
    fit3.SetParLimits(1,FitPeakPosX-1,FitPeakPosX+1)
    fit3.SetParLimits(2,FitSigma-1,FitSigma+1)
    fit3.SetParLimits(3,-1,1)
    fit3.SetParLimits(4,-2,2)
    fit3.SetParameters(FitHeight,FitPeakPosX,FitSigma,0,0) 
    HistClone.Fit(fit3,"Q","",peakPosX[i]-sigma,peakPosX[i]+sigma)
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
        canvas.Print(fitdir+filename+strPipeLab+str(i)+"_FitResults.png")

    #Height Finding
    fit4 = ROOT.TF1("pgaus4","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",Xmin,Xmax)
    fit4.SetParLimits(0,FitHeight-1,FitHeight+2)
    fit4.SetParLimits(1,FitPeakPosX-1,FitPeakPosX+1)
    fit4.SetParLimits(2,FitSigma-0.5,FitSigma+0.5)
    fit4.SetParLimits(3,-50,50)
    fit4.SetParLimits(4,-2,2)
    fit4.SetParameters(FitHeight,FitPeakPosX,FitSigma,0,0) 
    HistClone2.Fit(fit4,"Q","",peakPosX[i]-sigma,peakPosX[i]+sigma)
    FitHeightMeas = fit4.GetParameter(0)

    AcceptX = 2         #   Any peak that is found that is not +/- acceptace from
                        # the edges of the histogram or the peakX position 
                        # found by TSpectrum are discarded
    SizeConstant = 2.5  #   To be at 95% of the peaks area, one uses 4 
                        # (This seems to be about right for the flaw size 
                        # to implimented size ratio)

    TScale = 1
    if bolTempIsHot == False:
      TScale = 2

    if FitPeakPosX > peakPosX[i] -AcceptX and FitPeakPosX < peakPosX[i]+AcceptX and\
       FitPeakPosX>Xmin+AcceptX and FitPeakPosX< Xmax-AcceptX and\
       FitGoodness < 0.05*TScale:
      if abs(FitSigma*SizeConstant)<8 and abs(FitSigma*SizeConstant) > 1 and FitHeightMeas > 0.25*TScale: 
        MFlaws = np.append(MFlaws,[FitLevel,FitPeakPosX,FitSigma,FitHeightMeas,FitGoodness])
  
  #PLOTTING---------------------------------------------------------------
  fitline = ROOT.TF1("offset","[0]",Xmin,Xmax)
  Hist.Fit(fitline,"Q N","",Xmin,Xmax)
  Offset = fitline.GetParameter(0)

  canvas.Update()
  Hist.Draw()

  #Histogram Y Range
  if Offset > 22:
    Hist.SetMaximum(Offset +4)
    Hist.SetMinimum(Offset -4)
  elif Offset < -10:
    Hist.SetMaximum(Offset +4)
    Hist.SetMinimum(Offset -4)
  elif Offset < 10:
    Hist.SetMaximum(Offset +4)
    Hist.SetMinimum(Offset -4)
  else:
    Hist.SetMaximum(Offset+4)
    Hist.SetMinimum(Offset-4)
    
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
    canvas.Print(outdir+filename+strPipeLab+"Tpeaks.root")
    canvas.Print(outdir+filename+strPipeLab+"Tpeaks.png")

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
  DefectInfo = np.append(DefectInfo,FindPeaks(inputfile,outdir,fitdir,C1,0,bolVerb)) 
  DefectsTop = len(Boxes)
  for i in range(DefectsTop):
    Boxes[i].Draw()

  C1.cd(2)
  DefectInfo = np.append(DefectInfo,FindPeaks(inputfile,outdir,fitdir,C1,1,bolVerb))
  try:
    for i in range(len(Boxes)-DefectsTop):
      Boxes[i+DefectsTop].Draw()
  except:
    print("SOMETHING WENT WRONG")

  filename =MakeFileName(inputfile)

  C1.Print(outdir+filename+"-AllFlaws.png")
  C1.Print(outdir+filename+"-AllFlaws.root")
  PrintDefectInfo(DefectInfo,-1)
  return DefectInfo

#------------------------------------------------------------------------------
#Hot and Cold Defect Comparison------------------------------------------------
def HnCComp(filehot,filecold,outdir,fitdir,Canvas,bolVerb):
  """
  Takes two sets of data, both should be of the same stave in the same orientation,
  while one should be at a high temperature and the other at a low temperature. It
  will then combine the results from each
  """
  DefHot  = GetDefects(filehot,outdir,fitdir,Canvas,bolVerb)
  DefCold = GetDefects(filecold,outdir,fitdir,Canvas,bolVerb)


  Nthings = 7
  NHotDefects = len(DefHot)/Nthings
  NColdDefects = len(DefCold)/Nthings

  try:
    DefHot = DefHot.reshape(NHotDefects,Nthings)
  except:
    DefHot = []
  try:
    DefCold =DefCold.reshape(NColdDefects,Nthings)  
  except:
    DefCold = []
  #Convert defects from each image to each pipe (top or bottom)
  TDefs = []
  BDefs = []
  for i in range(NHotDefects):
    if DefHot[i][1]== 0:
      TDefs = np.append(TDefs,[0,0,DefHot[i][2],DefHot[i][3],DefHot[i][4],DefHot[i][5],DefHot[i][6]]) 
    else:
      BDefs = np.append(BDefs,[0,1,DefHot[i][2],DefHot[i][3],DefHot[i][4],DefHot[i][5],DefHot[i][6]])
  for i in range(NColdDefects):
    if DefCold[i][1]==0:
      TDefs = np.append(TDefs,[1,0,DefCold[i][2],DefCold[i][3],DefCold[i][4],DefCold[i][5],DefCold[i][6]])
    else:
      BDefs = np.append(BDefs,[1,1,DefCold[i][2],DefCold[i][3],DefCold[i][4],DefCold[i][5],DefCold[i][6]])

  NTDefs = len(TDefs)/Nthings
  NBDefs = len(BDefs)/Nthings

  print(TDefs)
  print(BDefs)

  try:
    TDefs = TDefs.reshape(NTDefs,Nthings)
    TDefs = TDefs[TDefs[:,2].argsort()]
  except:
    TDefs = []
  try:
    BDefs = BDefs.reshape(NBDefs,Nthings)
    BDefs = BDefs[BDefs[:,2].argsort()] 
  except:
    BDefs = []

  #Sorts the Defects by their position on the x axis

  Range = 1
  TDefsF = []
  TDefsO = [0,0,0,0,0,0,0]
  
  #Combines double counted defects to their average
  for i in range(NTDefs):
    if TDefsO[2] == 0 and i != NTDefs-1:
      pass
    elif TDefs[i][2] > TDefsO[2]-Range and TDefs[i][2] < TDefsO[2] + Range:
      TDefsF = np.append(TDefsF,[0.5*(TDefsO[0]+TDefs[i][0]),0,0.5*(TDefsO[2]+TDefs[i][2]),0.5*(TDefsO[3]+TDefs[i][3]),0.5*(TDefsO[4]+TDefs[i][4]),0.5*(TDefsO[5]+TDefs[i][5]),0.5*(TDefsO[6]+TDefs[i][6])])
      TDefsO[2] =0
      continue
    elif i == NTDefs-1:
      if TDefsO[2] != 0:
        TDefsF = np.append(TDefsF,TDefsO)
      TDefsF = np.append(TDefsF,TDefs[i])
    else:
      TDefsF = np.append(TDefsF,TDefsO)
    TDefsO = TDefs[i]
  BDefsF = []
  BDefsO = [0,0,0,0,0,0,0]
  for i in range(NBDefs):
    if BDefsO[2] == 0 and i != NBDefs-1:
      pass
    elif BDefs[i][2] > BDefsO[2]-Range and BDefs[i][2] < BDefsO[2] + Range:
      BDefsF = np.append(BDefsF,[0.5*(BDefsO[0]+BDefs[i][0]),1,0.5*(BDefsO[2]+BDefs[i][2]),0.5*(BDefsO[3]+BDefs[i][3]),0.5*(BDefsO[4]+BDefs[i][4]),0.5*(BDefsO[5]+BDefs[i][5]),0.5*(BDefsO[6]+BDefs[i][6])])
      BDefsO[2] =0
      continue
    elif i == NBDefs-1:
      if BDefsO[2] != 0:
        BDefsF = np.append(BDefsF,BDefsO)
      BDefsF = np.append(BDefsF,BDefs[i])
    else:
      BDefsF = np.append(BDefsF,BDefsO)
    BDefsO = BDefs[i]

  NTDefsF = len(TDefsF)/Nthings
  NBDefsF = len(BDefsF)/Nthings
 
  #Recombines all defects and their information together
  DefsF = np.append(TDefsF,BDefsF)
  print("COMBINED DEFECTS FOUND")
  PrintDefectInfo(DefsF)
  GetDefects2(filehot,filecold,outdir,fitdir,Canvas,bolVerb)

#------------------------------------------------------------------------------
def GetDefects2(inputfile1,inputfile2,outdir,fitdir,C1,bolVerb = 0):
  """
  This takes tow input files and gets all of the defects and prints all of the 
  relevant information. This is necessary for printing the boxes on a single canvas
  """
  DefectInfo = []
  C1.Clear()
  C1.Divide(1,4)
  C1.cd(1)
  global Boxes
  Boxes = []
  DefectInfo = np.append(DefectInfo,FindPeaks(inputfile1,outdir,fitdir,C1,0,bolVerb))
  C1.cd(2) 
  DefectInfo = np.append(DefectInfo,FindPeaks(inputfile2,outdir,fitdir,C1,0,bolVerb)) 
  C1.cd(3)
  DefectInfo = np.append(DefectInfo,FindPeaks(inputfile1,outdir,fitdir,C1,1,bolVerb))  
  C1.cd(4)
  DefectInfo = np.append(DefectInfo,FindPeaks(inputfile2,outdir,fitdir,C1,1,bolVerb))

  filename1 = MakeFileName(inputfile1)  
  filename2 = MakeFileName(inputfile2)

  C1.Print(outdir+filename1+filename2+"-AllFlaws.png")
  C1.Print(outdir+filename1+filename2+"-AllFlaws.root")

#------------------------------------------------------------------------------
def PrintDefectInfo(DefectInfo,intLine=-1):
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
      print " {0:2.1f}-{1:1}  {2:8.3f}  {3:8.3f}  {4:8.3f}  {5:8.5f}    {6:1}".format(DefectInfo[i*nThings],int(DefectInfo[i*nThings+1]),DefectInfo[i*nThings+2],DefectInfo[i*nThings+3],DefectInfo[i*nThings+4],DefectInfo[i*nThings+5],int(DefectInfo[i*nThings+6]))
  else:  
    print " {0:2.1f}-{1:1}  {2:8.3f}  {3:8.3f}  {4:8.3f}  {5:8.5f}    {6:1}".format(DefectInfo[intLine*nThings],int(DefectInfo[intLine*nThings+1]),DefectInfo[intLine*nThings+2],DefectInfo[intLine*nThings+3],DefectInfo[intLine*nThings+4],DefectInfo[intLine*nThings+5],int(DefectInfo[intLine*nThings+6]))
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
#Text Commands-----------------------------------------------------------------
def TextCommands(inputfile,outdir,fitdir,bolVerb):
  """
  This allows the user to give commands to replot items and keep the program and
  to edit what is plotted
  """
  print("\n-----Stave Defect Finder-----")
  print("h = print list of commands")
  print("q = quit program")
  print("f = shows all attached files")
  print("Defs = prints all defects from every file")
  print("Def(i) = prints ith inputfile")
  print("HnC(i,j) = prints hot and cold comparison between ith and jth inputfiles")
  print("TDiff(i,j,type,poff) = prints difference between ith and jth inputfiles type")
  print("                       where the type is given by temp, mean, width, chi2")
  print("                       and poff is an offset given to the first dataset")       
  print("RMS(i) = prints out a plot of the spectra fit with a line then subtracted")

  global bolQuit
  while bolQuit == False:
    Vin = raw_input("\nType the command you wish to do:\n")
    if Vin == 'q':
      bolQuit = True
    elif Vin == 'f':
      for i in range(len(inputfile)):
        print(str(i)+" "+inputfile[i])
    elif Vin == 'Defs':
      Cnew = ROOT.TCanvas("cnew","cnew",0,0,2000,600)
      Defs = []
      for i in inputfile:
        Defs = np.append(Defs,GetDefects(i,outdir,fitdir,Cnew,bolVerb))
        if bolVerb == True:
          DefectAnalysis(Def,outdir,Cnew)

    elif Vin == 'h':
      print("h = print list of commands")
      print("q = quit program")
      print("f = shows all attached files")
      print("Defs = prints all defects from every file")
      print("Def(i) = prints ith inputfile")
      print("HnC(i,j) = prints hot and cold comparison between ith and jth inputfiles")
      print("TDiff(i,j,type,poff) = prints difference between ith and jth inputfiles type")
      print("                       where the type is given by temp, mean, width, chi2")
      print("                       and poff is an offset given to the first dataset")       
      print("RMS(i) = prints out a plot of the spectra fit with a line then subtracted")

    elif 'Def(' in Vin:
      Vin = Vin.lstrip('Def(')
      Vin = Vin.rstrip(')')
      try:
        Vin = int(Vin)
        Cnew1 = ROOT.TCanvas("cnew1","cnew1",0,0,2000,600)
        Def = GetDefects(inputfile[Vin],outdir,fitdir,Cnew1,bolVerb)
      except:
        print("NOT Correct format") 
    elif 'HnC(' in Vin:
      Vin = Vin.lstrip('HnC(')
      Vin = Vin.rstrip(')')
      Vin = Vin.split(',')
      Cnew2 = ROOT.TCanvas("cnew2","cnew2",0,0,2000,1200)
      F = HnCComp(inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdir,fitdir,Cnew2,bolVerb)
      print("NOT Correct format")

    elif 'TDiff(' in Vin:
      Vin = Vin.lstrip('TDiff(')
      Vin = Vin.rstrip(')')
      Vin = Vin.split(',')
      ninputs = len(Vin)

      Ctop = ROOT.TCanvas("ctop","ctop",0,0,2000,490)
      Cbot = ROOT.TCanvas("cbot","cbot",0,510,2000,490)
      if ninputs == 2:
        try:
          F = PlotDiff(0,inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdir,Ctop,Cbot,'temperature;1')
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
          F = PlotDiff(0,inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdir,Ctop,Cbot,choice)
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
          F = PlotDiff(int(Vin[3]),inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdir,Ctop,Cbot,choice)
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
            
          F = PlotDiff(int(Vin[2]),inputfile[int(Vin[0])],inputfile[int(Vin[1])],outdir,Ctop,Cbot,options[i])          
      else:
        print("Not correct number of parameters")
    elif 'RMS(' in Vin:
      Vin = Vin.lstrip('RMS(')
      Vin = Vin.rstrip(')')
      try:
        Vin = int(Vin)
      except:
        continue
      CRMS = ROOT.TCanvas("crms","crms",0,0,2000,500)
      F = GetBothRMS(inputfile[Vin],outdir,CRMS)
      
        #print("Failed to Make RMS plot")
    else:
      print("Command not recognized")
#Find RMS---------------------------------------------------------------------
def findRMS(inputfile,intPipeNum,outdir):
  """
  Takes in a histogram, fits it with a strait line, it then subtracts each data
  point with the line. This is then plotted on a separate histogram.
  """
  #Get histogram info
  hist = GetHistogram(inputfile,intPipeNum)
  Info = GetHistInfo(hist)
  nbins = Info[0]
  Xmin = Info[1]
  Xmax = Info[2]

  #Fit histogram
  fit = ROOT.TF1("pol1","pol1",Xmin,Xmax)
  hist.Fit(fit,"Q","",Xmin,Xmax)
  offset =fit.GetParameter(0)
  slope =fit.GetParameter(1)
 
  #Create rms histogram 
  PipeName = "TopPipe"
  if intPipeNum >0:
    PipeName = "BotPipe"
  name = MakeFileName(inputfile)
  RMSHist = ROOT.TH1F(name+PipeName+"RMSHistogram",PipeName+" RMS Plot;Measured - Fit [#circC];Counts",50,-2,2)
  Xincrement = (Xmax - Xmin)/nbins
  for i in range(nbins):
    Val = hist.GetBinContent(i+1)
    LineVal = offset + slope*(Xmin + Xincrement*i)
    Diff = Val - LineVal
    RMSHist.Fill(Diff) 
    #print(Diff)

  return RMSHist

def GetBothRMS(inputfile,outdir,canvas):
  """
  Creates a plot that has both histograms on it
  """
  canvas.Divide(2,1)
  canvas.cd(1)
  h1 = findRMS(inputfile,0,outdir)
  canvas.Update()
  h1.Draw()
  canvas.cd(2)
  h2 = findRMS(inputfile,1,outdir)
  canvas.Update()
  h2.Draw()

  name = MakeFileName(inputfile)
  canvas.Print(outdir+name+"RMSPlot.png")

#------------------------------------------------------------------------------
#Main Routine------------------------------------------------------------------
def main():
  """
  The main routine
  """
  if sys.version_info[0] >= 3:
    print ("ERROR: PyROOT only works with Python 2.x")
    raise Exception(" Python Version too high. Use 2.x.")

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

  #MAKING THE OUTDIRECTORY
  outdir = "flawplots/"
  fitdir = outdir+"fits/"
  if not os.path.isdir( outdir ):
    os.mkdir( outdir )
  if not os.path.isdir( fitdir ):
    os.mkdir( fitdir )

  C1 = ROOT.TCanvas("c1","c1",0,0,2000,600)

  #First Run
  if len(inputfile) != 2:
    Def = GetDefects(inputfile[0],outdir,fitdir,C1,bolVerb)
    if bolVerb > 0:
      DefectAnalysis(Def,outdir,C1)
  else:
    Def = HnCComp(inputfile[0],inputfile[1],outdir,fitdir,C1,bolVerb)
   
  #Wait for commands
  global bolQuit 
  bolQuit = False
  while bolQuit == False:
    TextCommands(inputfile,outdir,fitdir,bolVerb)

if __name__ == "__main__":
  main()
