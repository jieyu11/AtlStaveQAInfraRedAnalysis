"""
This is the toolbox of functions for the defectFinder.py
"""
import sys          # system 
import os           # operating system
import ROOT         # ROOT from CERN
import numpy as np  # number in python


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
  strOrigName = Hist.GetName()
  strNamePrecurse = MakeFileName(inputfile)
  Hist.SetName(strNamePrecurse+'_'+strOrigName)
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
def InvertHistogram(HistToInvert):
  """
  Takes a histogram and inverts it around the y axis and returns the inverted hist
  """
  Info = GetHistInfo(HistToInvert)
  nbins = int(Info[0])
  Xmin = float(Info[1])
  Xmax = float(Info[2])
  HistInverted = ROOT.TH1F(Info[3]+"Inverted",Info[3]+"Inverted",nbins,Xmin,Xmax)
  for i in range(nbins):
    HistInverted.AddBinContent(i+1,HistToInvert.GetBinContent(i+1)*-1)
  return HistInverted

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
    return -999
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
def PlotDiffLine(offset,file1,file2,outdir,objCanvas,intPipeNum,strType,Scale = False):
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

  #Define Fits for Difference Plot
  nbins = Hist1.GetNbinsX()
  fit = ROOT.TF1( "FitLine", "pol0", 0., float(nbins) )
  Xmin= Hist1.GetXaxis().GetBinUpEdge(0)
  Xmax= Hist1.GetXaxis().GetBinUpEdge(nbins)

  Hist0Shifted = Hist0
  Hist1Shifted = Hist1

  #Rescaling the second plot
  if Scale == True:
    AmbTemp = 22
    Hist0Shifted.Fit(fit,"q nodraw")  
    Hist0Offset = fit.GetParameter(0)
    Hist1Shifted.Fit(fit,"q nodraw")
    Hist1Offset = fit.GetParameter(0)
    ScaleVal = Hist0Offset/Hist1Offset
    print(ScaleVal)
    for i in range(nbins):
      Hist1Shifted.SetBinContent(i+1,Hist1.GetBinContent(i+1)*ScaleVal)
  else:
    Hist1Shifted = Hist1

  #Fit the difference between the two plots
  HistDiff = Hist1Shifted - Hist0Shifted
  HistDiff.Fit(fit,"q")
  offset = fit.GetParameter(0)
  
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

  #Make Titles
  Hist0Title = MakeFileName(file1)
  Hist1Title = MakeFileName(file2)
  if Scale == True:
    Hist1Title = Hist1Title+"_scaled("+str(round(ScaleVal,3))+")"
    ScaleTitle = "SC"
  else:
    ScaleTitle = ""

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

  objCanvas.Print(outdir+filename+strPipeLab+strType[:-2]+ScaleTitle+"_dif.pdf")  
  objCanvas.Print(outdir+filename+strPipeLab+strType[:-2]+ScaleTitle+"_dif.root")
 

#------------------------------------------------------------------------------
def PlotDiff(offset,file1,file2,outdir,C1,C2,strType = "temperature;1"):
  """
  Plots the difference between two stave image temperatures along both the top and
  bottom cooling pipe.
  """

  C1.cd()
  C1.SetGrid()
  PlotDiffLine(offset,file1,file2,outdir,C1,0,strType)
  
  C2.cd() 
  C2.SetGrid()
  PlotDiffLine(offset,file1,file2,outdir,C2,1,strType) 

#------------------------------------------------------------------------------
def PlotDiffScale(offset,file1,file2,outdir,C1,C2,strType = "temperature;1"):
  """
  Plots the difference between two stave image temperatures along both the top and
  bottom cooling pipe.
  """

  C1.cd()
  C1.SetGrid()
  PlotDiffLine(offset,file1,file2,outdir,C1,0,strType,True)
  
  C2.cd()
  C2.SetGrid() 
  PlotDiffLine(offset,file1,file2,outdir,C2,1,strType,True) 

#-----------------------------------------------------------------------------
def GetBothRMS(inputfile,outdir,canvas):
  """
  Creates a plot that has both histograms on it
  """
  canvas.Divide(2,1)
  canvas.cd(1)
  HistTop = GetHistogram(inputfile,0)
  h1 = RMS(HistTop,outdir,canvas)
  canvas.Update()
  h1.Draw()
  canvas.cd(2)
  HistBot = GetHistogram(inputfile,1)
  h2 = RMS(HistBot,outdir,canvas)
  canvas.Update()
  h2.Draw()

  name = MakeFileName(inputfile)
  canvas.Print(outdir+name+"RMSPlot.pdf")
#Find RMS 2-------------------------------------------------------------------
def OneLineRMS(inputfile,outdir,canvas,bendlength = 10):
  """
  Stuff will go here
  """  
  canvas.Clear()
  #Create the line
  LineHist = OneLine(inputfile,outdir,canvas,bendlength)
  RMSHist = RMS(LineHist,outdir,canvas)
  canvas.cd()

  #Plot it
  RMSHist.Draw()
  canvas.Update()
  
  #Save it
  name = MakeFileName(inputfile)
  canvas.Print(outdir+name+"RMSPlotLine.pdf")


def RMS(hist,outdir,canvas):
  """
  Takes in a histogram, fits it with a strait line, it then subtracts each data
  point with the line. This is then plotted on a separate histogram.
  """
  #Get histogram info
  Info = GetHistInfo(hist)
  nbins = Info[0]
  Xmin = Info[1]
  Xmax = Info[2]
  name = Info[3]

  ROOT.gStyle.SetOptStat(1)
  ROOT.gStyle.SetOptTitle(1)

  #Fit histogram
  fit = ROOT.TF1("pol1","pol1",Xmin,Xmax)
  hist.Fit(fit,"Q","",Xmin,Xmax)
  offset =fit.GetParameter(0)
  slope =fit.GetParameter(1)
 
  #Create rms histogram  
  RMSHist = ROOT.TH1F(name+"RMSHistogram"," RMS Plot;Measured - Fit [#circC];Counts",50,-2,2)
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
  canvas.Print(outdir+name+"RMSPlot.pdf")


#-----------------------------------------------------------------------------
#One Line---------------------------------------------------------------------
def OneLine(inputfile,outdir,canvas,bendLength = 1):
  """
  Takes both lines of the pipe and stitches them together

  """
  ROOT.gStyle.SetOptStat(0)
  canvas.Clear()

  #Get histogram EOS info
  histEOS = GetHistogram(inputfile,0)
  InfoEOS = GetHistInfo(histEOS)
  nbinsEOS = InfoEOS[0]
  XminEOS = InfoEOS[1]
  XmaxEOS = InfoEOS[2]

  #Get histogram info
  hist1 = GetHistogram(inputfile,1)
  Info1 = GetHistInfo(hist1)
  nbins1 = Info1[0]
  Xmin1 = Info1[1]
  Xmax1 = Info1[2]

  #Finding Common Middle Line 
  LEnd = histEOS.GetBinContent(nbinsEOS)
  REnd = hist1.GetBinContent(nbins1)

  AvgOffset = 0.5*(LEnd+REnd)
  #print(AvgOffset)

  #Combining the two
  name = MakeFileName(inputfile)
  Conversion = (XmaxEOS-XminEOS)/nbinsEOS
  nbinsComb = nbinsEOS + nbins1 + int(bendLength/Conversion)
  Xmax = XmaxEOS +Xmax1-Xmin1+bendLength
  
  CombHist = ROOT.TH1F(name+"_BothPipes","Cooling Pipe Temp;Distance Along Pipe [cm];Temperature [#circC];",nbinsComb,XminEOS,Xmax)
  nbinsAfter = nbinsComb -nbins1
  for i in range(nbinsComb):
    if i < nbinsEOS:
      Val = histEOS.GetBinContent(i+1)
    elif i >= nbinsAfter:
      Val = hist1.GetBinContent(nbinsComb-i) 
    else: 
      Val = AvgOffset
    CombHist.SetBinContent(i+1,Val)

  #Create the plot
  #canvas.cd()
  #Pad0 = ROOT.TPad("grid","",0,0,1,1)
  #Pad0.Draw()
  #Pad0.cd()
  #Pad0.SetGrid()
  #CombHist.Draw()    
  #canvas.Update()
  #canvas.Print(outdir+name+"Line.root")
  #canvas.Print(outdir+name+"Line.pdf")

  return CombHist
  

#One Line Comp---------------------------------------------------------------------
def OneLineComp(inputfile0,inputfile1,outdir,canvas,bendLength = 1,Scale = False):
  """
  Takes two inputfiles and does the oneline function on them and then compares the two
  resulting histograms
  """
  ROOT.gStyle.SetOptStat(0)
  List = ROOT.TList()
  Hist0 = OneLine(inputfile0,outdir,canvas,bendLength)
  Hist1 = OneLine(inputfile1,outdir,canvas,bendLength)
  canvas.Clear()

  Hist0info = GetHistInfo(Hist0)
  Hist1info = GetHistInfo(Hist1)
  
  #Check to see if histograms are not on the same boundaries, fix them if they are.
  while Hist0info[0:2] != Hist1info[0:2]:
    Hists = CutHistToSameSize(Hist0,Hist1)
    Hist0 = Hists[0]
    Hist1 = Hists[1]
    Hist0info = GetHistInfo(Hist0)
    Hist1info = GetHistInfo(Hist1)

  #Define Fits for Difference Plot
  nbins = Hist1.GetNbinsX()
  fit = ROOT.TF1( "FitLine", "pol0", 0., float(nbins) )
  Xmin= Hist1.GetXaxis().GetBinUpEdge(0)
  Xmax= Hist1.GetXaxis().GetBinUpEdge(nbins)

  Hist0Shifted = Hist0
  Hist1Shifted = Hist1

  #Rescaling the second plot
  if Scale == True:
    Hist0Shifted.Fit(fit,"q nodraw")  
    Hist0Offset = fit.GetParameter(0)
    Hist1Shifted.Fit(fit,"q nodraw")
    Hist1Offset = fit.GetParameter(0)
    ScaleVal = Hist0Offset/Hist1Offset
    print(ScaleVal)
    for i in range(nbins):
      Hist1Shifted.SetBinContent(i+1,Hist1.GetBinContent(i+1)*ScaleVal)
  else:
    Hist1Shifted = Hist1

  #Fit the difference between the two plots
  HistDiff = Hist1Shifted - Hist0Shifted
  HistDiff.Fit(fit,"q nodraw")
  offset = fit.GetParameter(0)
  
  yMax = max(Hist0Shifted.GetMaximum(),Hist1Shifted.GetMaximum())
  yMin = min(Hist0Shifted.GetMinimum(),Hist1Shifted.GetMinimum())
 
  ySep = yMax-yMin

  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)

  #Pad Basics
  LS = 0.15
  TitleS = 0.15
  TitleO =0.7

  Pad1 = ROOT.TPad("pad1","Both Plots",0.0,0.3,1.0,1.0)
  Pad2 = ROOT.TPad("pad2","Difference Value",0.0,0.0,1.0,0.3)
  Pad1.Draw()
  Pad1.SetGrid()
  Pad2.Draw()
  Pad2.SetGrid()

  #Make Titles
  Hist0Title = MakeFileName(inputfile0)
  Hist1Title = MakeFileName(inputfile1)
  if Scale == True:
    Hist1Title = Hist1Title+"_scaled("+str(round(ScaleVal,3))+")"
    ScaleTitle = "SC"
  else:
    ScaleTitle = ""

  #Plot 1 Comparison
  Pad1.cd()
  Hist1Shifted.SetLineColor(1)
  Hist0Shifted.SetAxisRange(yMin-0.1*ySep,yMax+0.1*ySep,"Y")  
  Hist1Shifted.SetAxisRange(yMin-0.1*ySep,yMax+0.1*ySep,"Y") 
  Hist0Shifted.GetYaxis().SetLabelSize(LS*.3)
  Hist0Shifted.GetYaxis().SetTitleSize(TitleS*.3)
  Hist0Shifted.GetYaxis().SetTitleOffset(TitleO)
  Hist0Shifted.GetXaxis().SetLabelSize(LS*.3)
  Hist0Shifted.GetXaxis().SetTitleSize(TitleS*.3)
  Hist0Shifted.GetXaxis().SetTitleOffset(TitleO)

  Legend = ROOT.TLegend(0.1,0.7,0.3,0.9)
  Title = "Plot Comparison"
  Legend.SetHeader(Title,"C")
  Legend.AddEntry(Hist0Shifted,Hist0Title,"L")
  Legend.AddEntry(Hist1Shifted,Hist1Title,"L")

  List.AddLast(Hist0Shifted)
  List.AddLast(Hist1Shifted)
  List.AddLast(Legend)
  List.At(0).Draw("LC CM")
  List.At(1).Draw("SAME LC CM")
  List.At(2).Draw()

  #Plot 1 Delta
  Pad2.cd()
  HistDiff.SetDirectory(0) 
  HistDiff.GetYaxis().SetTitle("Temp Difference [C#circ]")
  HistDiff.GetYaxis().SetLabelSize(LS*.7)
  HistDiff.GetYaxis().SetTitleSize(LS*.7) 
  HistDiff.GetYaxis().SetTitleOffset(TitleO*0.5)
  HistDiff.GetXaxis().Delete()
  HistDiff.SetNdivisions(10,"Y")
  HistDiff.GetXaxis().SetLabelOffset(10)

  List.AddLast(HistDiff)
  List.At(3).Draw("LC CM")
 
  #PRINTING THE PLOTS    
  filename = Hist0Title+Hist1Title 
  canvas.Print(outdir+filename+ScaleTitle+"_LineDif.pdf")  
  canvas.Print(outdir+filename+ScaleTitle+"_LineDif.root") 
#---------------------------------------------------------------------
def FindTempHist(inputfile,outdir,canvas):
  """
  Plots all temperature measurements of the profile as a histogram
  """
 
  #Get histogram EOS info
  histEOS = GetHistogram(inputfile,0)
  InfoEOS = GetHistInfo(histEOS)
  nbinsEOS = InfoEOS[0]
  XminEOS = InfoEOS[1]
  XmaxEOS = InfoEOS[2]
  YmeanEOS = 0.5*(histEOS.GetMaximum()+histEOS.GetMinimum())

  #Get histogram info
  hist1 = GetHistogram(inputfile,1)
  Info1 = GetHistInfo(hist1)
  nbins1 = Info1[0]
  Xmin1 = Info1[1]
  Xmax1 = Info1[2]
  Ymean1 = 0.5*(hist1.GetMaximum()+hist1.GetMinimum())

  #Make new Histogram
  Ymean = 0.5*(YmeanEOS+Ymean1)
  Ysep = 10
  Ymin = Ymean -0.5*Ysep
  Ymax = Ymean +0.5*Ysep  

  HistNet = ROOT.TH1F("HistogramSum","HistSum",100,Ymin,Ymax)
  for i in range(nbinsEOS):
    HistNet.Fill(histEOS.GetBinContent(i+1))
    HistNet.Fill(histEOS.GetBinContent(i+1))
  
  #Print out the mean, std dev, and plot
  MeanVal = HistNet.GetMean()
  MeanErr = HistNet.GetStdDev()

  print("Name  : "+inputfile)
  print("Mean  : "+str(MeanVal))
  print("StdDev: "+str(MeanErr))

  nametag = MakeFileName(inputfile)
  ROOT.gStyle.SetOptStat(1)
  canvas.cd()
  canvas.Clear()
  HistNet.Draw()
  canvas.Print(outdir+nametag+'TempHist.pdf')

  return [MeanVal,MeanErr]

def FindAvgTemps(inputfiles,outdir,canvas):
  """
  Prints out a cvs file of the data as a function of a set time
  """
  outputFile = open('AvgTempData.csv','w')
  StartLine = 'Name,Mean,StdDev\n'
  outputFile.write(StartLine)
  outputFile.close()
  outputFile = open('AvgTempData.csv','a')
  counter = 0
  for fyle in inputfiles:
    Data = FindTempHist(fyle,outdir[counter],canvas)
    strData = str(Data[0])+','+str(Data[1])
    outline = fyle+','+strData+'\n'
    outputFile.write(outline)
  outputFile.close
 
#---------------------------------------------------------------------
def OneLineMulti(inputfiles,outdir,canvas,bendLength = 12):
  """
  Plots all attached input files on one plot

  WORKS WELL (March 5, 2019)
  """

  Hists = ROOT.TList() #Place to store all of the plots
  YMax =-999
  YMin = 999
  canvas.Clear()
  canvas.cd()
  canvas.SetGrid(1,1)
  Legend = ROOT.TLegend(0.1,0.7,0.3,0.9)
  Legend.SetHeader("Temperature Comparison","C")
  Hists.AddLast(Legend)
  name = ""

  for i in range(len(inputfiles)):
    HistN = OneLine(inputfiles[i],outdir[i],canvas,bendLength)
    YMax = max(YMax,HistN.GetMaximum())
    YMin = min(YMin,HistN.GetMinimum())
    name += HistN.GetName()
    HistN.SetAxisRange(YMin,YMax,"Y")
    HistN.SetLineColor(i+1)
    Hists.AddLast(HistN) 
    Legend.AddEntry(HistN,HistN.GetName())
  for i in range(1,len(inputfiles)+1):
    if i == 1: #Reset Base Drawing
      Hists.At(i).Draw("LC CM")
      Hists.At(i).SetAxisRange(YMin,YMax,"Y")
    else:
      Hists.At(i).Draw("Same LC CM")
      Hists.At(i).SetAxisRange(YMin,YMax,"Y")
  Legend.Draw()
  if len(name) > 30:
    listname = name.split("Rec-")
    NewName = ""
    for i in listname:
      NewName += i[0:6]+"_"
    name = NewName
  canvas.Update()
  canvas.Print(outdir[0]+name+"CombinedTempProfile.pdf")

  return
  #Makes frames for a video!
  canvas2 = ROOT.TCanvas("vidplot")

  j=-4.48
  for i in range(1,len(inputfiles)+1):
    canvas2.Clear()
    canvas2.cd()
    TimeStamp = ROOT.TPaveText(0.7,0.8,.9,.9,"NDC")
    j+=1

    TimeStamp.AddText("Time = "+str(j)+" min")
    Hists.At(i).SetLineColor(1)
    Hists.At(i).Draw("LC CM")
    Hists.At(i).SetAxisRange(YMin,YMax,"Y")
    TimeStamp.Draw()
    canvas2.Print(outdir[0]+Hists.At(i).GetName()+"frame.pdf")
 
