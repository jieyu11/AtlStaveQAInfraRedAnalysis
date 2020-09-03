"""
This is a portion of the defectFinder.py code that finds defects

"""

import sys          # system 
import os           # operating system
import ROOT         # ROOT from CERN
import numpy as np  # number in pyth
from defectFinderToolBox import *


#------------------------------------------------------------------------------
#FindPeaks---------------------------------------------------------------------

def FindPeaks(Hist,outdir,fitdir,canvas,bolVerb = 0): 
  """
  Uses the information along the cooling pipe histogram to find the flaws along the pipe.
  It has three main parts.

  1. Histogram conversion:
    - Uses input histogram
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
    Each peak that passes the cut is then saved as a defect. The flaw information
    is then returned by the function.
  """
  
  Info = GetHistInfo(Hist)
  nBins = Info[0] 
  Xmin = Info[1] 
  Xmax =  Info[2] 
  filename = Info[3]

  npeaks = 20
  objSpectrum = ROOT.TSpectrum(2*npeaks)
  HistClone =  Hist.Clone()
 
  bolTempIsHot = TempIsHot(Hist)
  if bolTempIsHot == True:
    for i in range(nBins):
      HistClone.SetBinContent(i+1,-1*Hist.GetBinContent(i+1)) 
  HistClone2 = HistClone.Clone() 

  #FINDING FLAWS---------------------------------------------------------- 
  HistClone.Draw()
  if bolVerb > 0:
    canvas.Print(fitdir+filename+"InvertPlot.pdf")
 


  Hist_SM_IV = HistClone.Clone()#BandPassFFT(HistClone,2,200) 

  Hist_SM_IV.Draw() 
  if bolVerb > 0:
    canvas.Print(fitdir+filename+"SmoothedInvertPlot.pdf")


  HistWhatsRemoved = HistClone2 - Hist_SM_IV
  HistWhatsRemoved.Draw() 
  if bolVerb > 0:
    canvas.Print(fitdir+filename+"WhatsRemovedPlot.pdf")
 
  HistBackground = objSpectrum.Background(Hist_SM_IV,17,"")

  HistPeaks = Hist_SM_IV - HistBackground

  sigma = 8
  threshold = 0.05
  nfound = objSpectrum.Search(HistPeaks,sigma,"noMarkov",threshold) #This removes the background again... I have left it in
  peakPosX = objSpectrum.GetPositionX()

  Hist_SM_IV.Draw()
  HistBackground.Draw("same")
  if bolVerb > 0:
    canvas.Print(fitdir+filename+"BackgroundInvertPlot.pdf")
 
  PeakPosArray = np.zeros(nfound) 
  for i in range(nfound):
    PeakPosArray[i] = peakPosX[i]  
  #print (str(PeakPosArray))
  
  MFlaws = []
  for i in range(nfound):
    #Find the Flaws Method 1---------- #Find Peak Shape 
    fit = ROOT.TF1("pgaus1","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",Xmin,Xmax)
    fit.SetParLimits(0,0,2)
    fit.SetParLimits(1,peakPosX[i]-sigma/2,peakPosX[i]+sigma/2)
    fit.SetParLimits(2,0.25,10)
    fit.SetParameters(1,peakPosX[i],1,0) 
    HistPeaks.Fit(fit,"Q","",peakPosX[i]-sigma/2,peakPosX[i]+sigma/2)
    #Get the fit info
    FitPeakPosX = fit.GetParameter(1)
    FitSigma = abs(fit.GetParameter(2))
    FitChiSq = fit.GetChisquare()
    FitDegFr = fit.GetNDF()
    FitHeight = abs(fit.GetParameter(0)) 
    FitGoodness = abs(FitChiSq/FitDegFr)
    FitLevel = 0
    if bolVerb > 0:
      canvas.Print(fitdir+filename+str(i)+"_FitResults0.pdf")

    #Find peak offset
    fitOff = ROOT.TF1("line","[0]",Xmin,Xmax)
    HistClone2.Fit(fitOff,"Q","")
    FitOffset = fitOff.GetParameter(0)

    #Find the Flaws Method 2---------- #Fit with constant background
    fit2 = ROOT.TF1("pgaus2","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]",Xmin,Xmax)
    fit2.SetParLimits(0,FitHeight,FitHeight*5)
    fit2.SetParLimits(1,peakPosX[i]-5,peakPosX[i]+5)
    fit2.SetParLimits(2,FitSigma-FitSigma*0.1,FitSigma+FitSigma*0.1)
    fit2.SetParLimits(3,FitOffset-5,FitOffset+5)
    fit2.SetParameters(FitHeight,FitPeakPosX,FitSigma,FitOffset) 
    HistClone2.Fit(fit2,"Q","",peakPosX[i]-4*FitSigma,peakPosX[i]+4*FitSigma)
    #Get the fit info
    Fit2PeakPosX = fit2.GetParameter(1)
    Fit2Sigma = abs(fit2.GetParameter(2))
    Fit2ChiSq = fit2.GetChisquare()
    Fit2DegFr = fit2.GetNDF()
    Fit2Height = abs(fit2.GetParameter(0)) 
    Fit2Goodness = abs(Fit2ChiSq/Fit2DegFr)

    #This fit should be worse... but it requires the starting conditions from the first fit
    #if Fit2Goodness < FitGoodness:
    FitPeakPosX = Fit2PeakPosX
    FitSigma = Fit2Sigma
    FitHeight = Fit2Height
    FitGoodness = Fit2Goodness
    FitLevel = 1
    if bolVerb > 0:
      canvas.Print(fitdir+filename+str(i)+"_FitResults1.pdf")

    #Find Flaws Method 3---------- #Fit with linear background
    fit3 = ROOT.TF1("pgaus3","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",Xmin,Xmax)
    fit3.SetParLimits(0,FitHeight-0.1,FitHeight+0.1)
    fit3.SetParLimits(1,FitPeakPosX-1,FitPeakPosX+1)
    fit3.SetParLimits(2,FitSigma-1,FitSigma+1)
    fit3.SetParLimits(3,FitOffset-10,FitOffset+10)
    fit3.SetParLimits(4,-2,2)
    fit3.SetParameters(FitHeight,FitPeakPosX,FitSigma,FitOffset,0) 
    HistClone2.Fit(fit3,"Q","",peakPosX[i]-4*FitSigma,peakPosX[i]+4*FitSigma)
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
      canvas.Print(fitdir+filename+str(i)+"_FitResults2.pdf")

    #Find Flaws Method 4---------- #Fit with parabolic background
    fit4 = ROOT.TF1("pgaus4","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x+[5]*x*x",Xmin,Xmax)
    fit4.SetParLimits(0,FitHeight-1,FitHeight+2)
    fit4.SetParLimits(1,FitPeakPosX-1,FitPeakPosX+1)
    fit4.SetParLimits(2,FitSigma-0.5,FitSigma+0.5)
    fit4.SetParLimits(3,FitOffset-10,FitOffset+10)
    fit4.SetParLimits(4,fit3.GetParameter(4)-2,fit3.GetParameter(4)+2)
    fit4.SetParLimits(5,-5,5)
    fit4.SetParameters(FitHeight,FitPeakPosX,FitSigma,FitOffset,fit3.GetParameter(4),0) 
    HistClone2.Fit(fit4,"Q","",peakPosX[i]-4*FitSigma,peakPosX[i]+4*FitSigma)
    try:
      Fit4Goodness = abs(fit4.GetChisquare()/fit4.GetNDF())
    except:
      Fit4Goodness = 999
    if Fit4Goodness < FitGoodness:
      FitPeakPosX = fit4.GetParameter(1)
      FitSigma = abs(fit4.GetParameter(2))
      FitHeight = abs(fit4.GetParameter(0))
      FitGoodness = Fit4Goodness
      FitLevel = 3
    if bolVerb > 0:
      canvas.Print(fitdir+filename+str(i)+"_FitResults3.pdf")

    FitHeightMeas = FitHeight

    AcceptX = 2         #   Any peak that is found that is not +/- acceptace from
                        # the edges of the histogram or the peakX position 
                        # found by TSpectrum are discarded
    SizeConstant = 2.5  #   To be at 95% of the peaks area, one uses 4 
                        # (This seems to be about right for the flaw size 
                        # to implimented size ratio)

    TScale = 1
    if bolTempIsHot == False:
      TScale = 2

    if FitPeakPosX > peakPosX[i] -AcceptX and\
       FitPeakPosX < peakPosX[i]+AcceptX and\
       FitPeakPosX>Xmin+AcceptX and\
       FitPeakPosX< Xmax-AcceptX and\
       FitGoodness < 0.05*TScale:
      if abs(FitSigma*SizeConstant)<8 and abs(FitSigma*SizeConstant) > 1 and FitHeightMeas > 0.25*TScale: 
        MFlaws = np.append(MFlaws,[FitLevel,FitPeakPosX,FitSigma,FitHeightMeas,FitGoodness])
  
  #PLOTTING---------------------------------------------------------------
  fitline = ROOT.TF1("offset","[0]",Xmin,Xmax)
  Hist.Fit(fitline,"Q N","",Xmin,Xmax)
  Offset = fitline.GetParameter(0)

  canvas.Update()
  Hist.Draw()

  LS = 0.15
  TitleO=0.7
  Hist.GetYaxis().SetLabelSize(LS*.7)
  Hist.GetYaxis().SetTitleSize(LS*.7) 
  Hist.GetYaxis().SetTitleOffset(TitleO*0.5)
  Hist.GetXaxis().SetLabelSize(LS*.4)
  Hist.GetXaxis().SetTitleSize(LS*.4)
  Hist.GetXaxis().SetTitleOffset(TitleO*1.2)
 
  #Histogram Y Range
  """
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
  """
  YMax = Hist.GetMaximum()
  YMin = Hist.GetMinimum()
  YSep = abs(YMax - YMin)
  PPer = 0.05
  Hist.SetMaximum(YMax+YSep*PPer)
  Hist.SetMinimum(YMin -YSep*PPer)  
  nMflaws = len(MFlaws)/5

  if bolVerb > 0: 
    canvas.Print(outdir+filename+"Tpeaks.root")
    canvas.Print(outdir+filename+"Tpeaks.pdf")

  #PREPPINGDEFECTINFO-----------------------------------------------------
  if 'top' in filename:
    intPipeNum = 0
  elif 'bot' in filename:
    intPipeNum = 1
  else:
    intPipeNum = -1

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
def GetDefectBoxes(DefectInfo,Ymin,Ymax,bolTempHot = True):
  """
  Creates a list of TBox objects for each defect
  """
 
  #Ymax = Hist.GetMaximum()
  #Ymin = Hist.GetMinimum()
  Boxes = []
  TScale = 1
  if bolTempHot == False:
    TScale = 2.
 #Major Flaw Boxes
  for i in range(len(DefectInfo)/7):
    Center = DefectInfo[i*7+2]
    Xmin = Center - 0.5*DefectInfo[i*7+3]
    Xmax = Center + 0.5*DefectInfo[i*7+3] 
    Boxes = np.append(Boxes,ROOT.TBox(Xmin,Ymin,Xmax,Ymax))
    FitHeight = DefectInfo[i*7+4]

    trans = 0.5
    if FitHeight > 1.6*TScale:
      Boxes[i].SetFillColorAlpha(2,trans)
    else:
      Boxes[i].SetFillColorAlpha(4,trans)
 
  return Boxes

#------------------------------------------------------------------------------
def GetDefects(inputfile,outdir,fitdir,C1,bolVerb = 0):
  """
  This takes an input file and gets all of the defects and prints all of them on 
  a single plot that is broken into the top(EOS) and bottom pipes thermal profile.
  It then returns the Defect Info 
  """
  #Find the Defects
  HistTop = GetHistogram(inputfile,0)
  DefectInfoTop = FindPeaks(HistTop,outdir,fitdir,C1,bolVerb)
  HistBot = GetHistogram(inputfile,1)
  DefectInfoBot = FindPeaks(HistBot,outdir,fitdir,C1,bolVerb)
  C1.Clear()
  C1.Update()

  Ymax = max(HistTop.GetMaximum(),HistBot.GetMaximum())
  Ymin = min(HistTop.GetMinimum(),HistBot.GetMinimum())

  bolTempHot = TempIsHot(HistTop)
  if bolTempHot == True:
    Scale = 5
  elif bolTempHot == False:
    Scale = 10
  else:
    Scale = 15

  YCent = (Ymax+Ymin)/2
  Ymax = YCent + Scale/2
  Ymin = YCent - Scale/2 

  #Plot the Histograms and Defects
  C1.Divide(1,2)
  C1.SetGrid()
  C1.cd(1)

  HistTop.Draw()
  HistTop.SetAxisRange(Ymin,Ymax,"Y")
  BoxesTop = GetDefectBoxes(DefectInfoTop,Ymin,Ymax,bolTempHot)
  for Box in BoxesTop:
    Box.Draw() 
  #HistTop.SetDirectory(0)
  C1.cd(2)
  HistBot.Draw()
  HistBot.SetAxisRange(Ymin,Ymax,"Y")
  BoxesBot = GetDefectBoxes(DefectInfoBot,Ymin,Ymax,bolTempHot)
  for Box in BoxesBot:
    Box.Draw()
  #HistBot.SetDirectory(0)
  C1.Update()
  DefectInfo = np.append(DefectInfoTop,DefectInfoBot)

  #Save the output as a pdf and root file
  filename = MakeFileName(inputfile)
  C1.Print(outdir+filename+"-AllFlaws.pdf")
  C1.Print(outdir+filename+"-AllFlaws.root")

  #Output the defect information
  PrintDefectInfo(DefectInfo,-1,bolTempHot)
  
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
  PlotDefectsHnC(filehot,filecold,TDefsF,BDefsF,outdir,fitdir,Canvas,bolVerb)

#------------------------------------------------------------------------------
def PlotDefectsHnC(inputfile1,inputfile2,TDefects,BDefects,outdir,fitdir,C1,bolVerb = 0):
  """
  This takes two input files and gets all of the defects and prints all of the 
  relevant information. This is necessary for printing the boxes on a single canvas
  """
  HistTop1 = GetHistogram(inputfile1,0)
  HistBot1 = GetHistogram(inputfile1,1)
  HistTop2 = GetHistogram(inputfile2,0)
  HistBot2 = GetHistogram(inputfile2,1)

  #Invert Cold spectra
  bolTempHot1 = TempIsHot(HistTop1)
  bolTempHot2 = TempIsHot(HistTop2)

  '''
  if bolTempHot1 == False:
    HistTop1 = InvertHistogram(HistTop1)
    HistBot1 = InvertHistogram(HistBot1)
  if bolTempHot2 == False:
    HistTop2 = InvertHistogram(HistTop2)
    HistBot2 = InvertHistogram(HistBot2)
  '''

  #Create Pads!
  C1.Clear()
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)

  PadT = ROOT.TPad("toppad","EOSPipe",0.0,0.5,1.0,1.0)
  PadB = ROOT.TPad("botpad","OutPipe",0.0,0.0,1.0,0.5)
  PadT.Draw()
  PadT.SetGrid()
  PadB.Draw()
  PadB.SetGrid()

  #Plotting Ranges
  Max = max(HistTop1.GetMaximum(),HistTop2.GetMaximum(),HistBot1.GetMaximum(),HistBot2.GetMaximum())
  Min = min(HistTop1.GetMinimum(),HistTop2.GetMinimum(),HistBot1.GetMinimum(),HistBot2.GetMinimum())
  YSep = Max-Min

  LS =0.15
  TitleS =0.15
  TitleO =0.5
  OffsetP = 0.5

  #Top Pad
  PadT.cd()
  HistTop1.Draw("LC CM")
  HistTop2.Draw("SAME LC CM")
  HistTop1.SetAxisRange(Min-0.1*YSep,Max+OffsetP*YSep,"Y")
  HistTop2.SetAxisRange(Min-0.1*YSep,Max+OffsetP*YSep,"Y")
  HistTop1.GetYaxis().SetLabelSize(LS*.3)
  HistTop1.GetYaxis().SetTitleSize(LS*.3)
  HistTop1.GetYaxis().SetTitleOffset(TitleO)
  HistTop1.GetYaxis().SetTitle("Temperature [#circC]")
  HistTop1.GetXaxis().SetLabelSize(LS*.3)
  HistTop1.GetXaxis().SetTitleSize(LS*.3)
  HistTop1.GetXaxis().SetTitleOffset(TitleO*1.5)
  HistTop1.GetXaxis().SetTitle("Stave Length [cm]")

  TLegend = ROOT.TLegend(0.1,0.7,0.3,0.9)
  TLegend.SetHeader("End of Stave Pipe Comparison","C")
  TLegend.AddEntry(HistTop1,HistTop1.GetName(),"L")
  TLegend.AddEntry(HistTop2,HistTop2.GetName(),"L")
  TLegend.Draw()

  TDefBoxes = GetDefectBoxes(TDefects,Min,Max,bolTempHot=True)
  for Box in TDefBoxes:
    Box.Draw()

  #Bottom Pad
  PadB.cd()
  HistBot1.Draw("LC CM")
  HistBot2.Draw("SAME LC CM")
  HistBot1.SetAxisRange(Min-0.1*YSep,Max+OffsetP*YSep,"Y")
  HistBot2.SetAxisRange(Min-0.1*YSep,Max+OffsetP*YSep,"Y")
  HistBot1.GetYaxis().SetLabelSize(LS*.3)
  HistBot1.GetYaxis().SetTitleSize(LS*.3)
  HistBot1.GetYaxis().SetTitleOffset(TitleO)
  HistBot1.GetYaxis().SetTitle("Temperature [#circC]")
  HistBot1.GetXaxis().SetLabelSize(LS*.3)
  HistBot1.GetXaxis().SetTitleSize(LS*.3)
  HistBot1.GetXaxis().SetTitleOffset(TitleO*1.5)
  HistBot1.GetXaxis().SetTitle("Stave Length [cm]")

  BLegend = ROOT.TLegend(0.1,0.7,0.3,0.9)
  BLegend.SetHeader("Return Pipe Comparison","C")
  BLegend.AddEntry(HistBot1,HistBot1.GetName(),"L")
  BLegend.AddEntry(HistBot2,HistBot2.GetName(),"L")
  BLegend.Draw()


  BDefBoxes = GetDefectBoxes(BDefects,Min,Max,bolTempHot=True)
  for Box in BDefBoxes:
    Box.Draw()
  C1.Update()

  #Save the plots
  filename1 = MakeFileName(inputfile1)  
  filename2 = MakeFileName(inputfile2)

  C1.Print(outdir+filename1+filename2+"-AllFlaws.pdf")
  C1.Print(outdir+filename1+filename2+"-AllFlaws.root")

#------------------------------------------------------------------------------
def SigDefect(peakHeight,bolTempHot):
  """
  This prints out a FLAWED Stave statement during defect info printing
  """

  TScale = 1.
  if bolTempHot == False:
    TScale = 2.
  if peakHeight >= 1.6*TScale:
    print("        !!!! SIGNIFICANT DEFECT FOUND !!!!")

 
#------------------------------------------------------------------------------
def PrintDefectInfo(DefectInfo,intLine=-1,bolTempHot=True):
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
      print("array index: {}".format(int(i*nThings+4)))
      print(DefectInfo)
      SigDefect(DefectInfo[int(i*nThings+4)],bolTempHot)
      print(" {0:2.1f}-{1:1}  {2:8.3f}  {3:8.3f}  {4:8.3f}  {5:8.5f}    {6:1}".format(DefectInfo[i*nThings],int(DefectInfo[i*nThings+1]),DefectInfo[i*nThings+2],DefectInfo[i*nThings+3],DefectInfo[i*nThings+4],DefectInfo[i*nThings+5],int(DefectInfo[i*nThings+6])))
      SigDefect(DefectInfo[i*nThings+4],bolTempHot)
  else: 
    SigDefect(DefectInfo[intLine*nThings+4],bolTempHot)
    print(" {0:2.1f}-{1:1}  {2:8.3f}  {3:8.3f}  {4:8.3f}  {5:8.5f}    {6:1}".format(DefectInfo[intLine*nThings],int(DefectInfo[intLine*nThings+1]),DefectInfo[intLine*nThings+2],DefectInfo[intLine*nThings+3],DefectInfo[intLine*nThings+4],DefectInfo[intLine*nThings+5],int(DefectInfo[intLine*nThings+6])))
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
    Canvas.Print(outdir+Name[i]+"_Hist.pdf")

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
  
  Canvas.Print(outdir+"WidthVHeightGraph.pdf")
  Canvas.Print(outdir+"WidthVHeightGraph.root")

#------------------------------------------------------------------------------
def GetOneLineDefects(inputfile,outdir,fitdir,canvas,bolVerb = 0):
  """
  This takes an input file and gets all of the defects and prints all of them on 
  a single plot that is broken into the top(EOS) and bottom pipes thermal profile.
  It then returns the Defect Info 
  """
  #Find the Defects
  bendlength = 12 #Not certain if this is what we want...
  OneLineHist = OneLine(inputfile,outdir,canvas,bendlength)
  DefectInfo = FindPeaks(OneLineHist,outdir,fitdir,canvas,bolVerb)
  
  canvas.Clear()
  canvas.Update()

  bolTempHot = TempIsHot(OneLineHist)
    
  AYMax = OneLineHist.GetMaximum()
  AYMin = OneLineHist.GetMinimum()
  if bolTempHot == True:
    Scale = 5
  elif bolTempHot == False:
    Scale = 10
  
  YCent = (AYMax + AYMin)/2
  YMax = YCent + Scale/2
  YMin = YCent - Scale/2

  #Plot the Histograms and Defects
  canvas.SetGrid()
  canvas.cd()

  OneLineHist.Draw()

  LS = 0.1
  TitleO = 1
  OneLineHist.SetAxisRange(YMin,YMax,"Y")
  OneLineHist.GetYaxis().SetLabelSize(LS*0.3)
  OneLineHist.GetYaxis().SetTitleSize(LS*.3)
  OneLineHist.GetYaxis().SetTitleOffset(TitleO)
  OneLineHist.GetXaxis().SetLabelSize(LS*0.3)
  OneLineHist.GetXaxis().SetTitleSize(LS*.3)
  OneLineHist.GetXaxis().SetTitleOffset(TitleO)


  Boxes = GetDefectBoxes(DefectInfo,YMin,YMax,bolTempHot)
  for Box in Boxes:
    Box.Draw() 
  canvas.Update()

  #Save the output as a pdf and root file
  filename = MakeFileName(inputfile)
  canvas.Print(outdir+filename+"-AllFlawsLine.pdf")
  canvas.Print(outdir+filename+"-AllFlawsLine.root")
 
  return DefectInfo


