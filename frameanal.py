#!/usr/bin/python

"""
@run:
  ./frameanal.py Test.seq

@brief:
  The code reads the root file produced for a IR camera temperature image, find
  the cooling pipe temperature profile which is used to detect the potential flaws
  on the stave.
  
@functions (class FrameAnalysis):
  __init__ (roo_name, cfg_name = "config_frame", fig_outdir = "plot") :
    Constructor providing the input root file name, roo_name; 
    config file name, cfg_name and output figure folder name, fig_outdir.
    In the config file, set the parameters:
    - Set the stave starting (X0, Y0) and ending (X1, Y1) pixel indexes, and 
      similarly for pipe region. 
    - Set the CM per pixel value, which depends on the IR camera and its lens. 
    - Set StaveSideL = 0 or 1 to indicate if the initial stave is on L (1) or 
      J (0) side. In the case of L side, it is flipped to look like J side.
    - Set LiquidTLow = 0 or 1 to indicate if the liquid passing the cooling 
      pipe is below (1) or above (0) room temperature.
    - Set the minimum and maximum temperature for plotting the raw frame, 
      stave region and pipe region.

  draw_frames():
    Draw the 2D temperature profile for the raw image, stave region and pipe region

  fit_hist(h1):
    Fit a gaussian to the input 1D histogram, h1. The X axis of h1 is the pixel number 
      from Y axis of the pipe region 2D image at a each X axis point. The Y axis of h1
      is the temperature values at the point.
    The fit returns a tuple: (temperature, mean position, width, chi2, ndf)
    The returned temperature is the measurement of the cooling pipe temperature at that
      point.

  find_pipes():
    Find the cooling pipes on the pipe region image. Make plots of the temperature
      distribution along the cooling pipe, together with the fitted result of the mean
      position, peak width, chi2 and NDF.

@notes:
  This code assumes using python 2.7.10
  It doesn't work for python 3.x, because ROOT support is not guaranteed.

@reference:
  TTree in python: 
    https://www-zeuthen.desy.de/~middell/public/pyroot/pyroot.html
    https://root.cern.ch/how/how-write-ttree-python

@email: jie.yu@cern.ch
"""

import sys   # system
import os    # operating system
import numpy # number in python
import math  # math
import ROOT  # ROOT from CERN

import configFinder as cf

class FrameAnalysis:
  """
    To find out the possible flaws on a stave, one needs to check the non-uniformity on the cooling pipe temperature
    profile.

    A stave recorded by IR camera is characterized with its frame diameter from bottom left corner to top right corner.
    It should look like:
    
      A raw "J" side IR camera 'photo'
   | -----------------------------------| 
   |                                    | 
   |                         Stave      |
   |     ---- .............. (X1, Y1)   |
   |     |  |                           |
   |     |   --------------- Pipe       |
   |     |-----------------| (X1, Y1)   |
   | Pipe|                 |            |
   | (X0,|-----------------|            |
   |  Y0)-------------------            |
   |     Stave                          |
   |     (X0, Y0)                       |
   |                                    |
   | -----------------------------------|  
   A raw "L" side image is flipped to look like "J" side above.

   A stave contains the whole stave including the end of stave card, indicated as Stave (X0, Y0) --> (X1, Y0)
   A pipe part of the stave contains the pipe region, indicated as Pipe (X0, Y0) --> (X1, Y0)
   The cooling pipe temperature is extracted from the Pipe frame.
    
  """
  _parameters = { 
    "StavePixelX0": 0, # stave position in pixel index
    "StavePixelY0": 0,
    "StavePixelX1": 0, 
    "StavePixelY1": 0,
    "PipePixelX0":  0, # pipe position in pixel index
    "PipePixelY0":  0,
    "PipePixelX1":  0, 
    "PipePixelY1":  0,
    "CMperPixel":   0, # cm per pixel, specific for camera
    "StaveSideL":   0, # 0: "J", 1: "L". Sides "J" or "L" based on the rotation of the end of stave card
    "LiquidTLow":   0, # 0: higher than room temperature, 1: lower than room temperature
    "FrameTmax":    -999., # for plotting, set the maximum and minimum Temperature
    "FrameTmin":     999., # use arbitrary value if these parameters are not set
    "StaveTmax":    -999., # Frame = raw frame, Stave = stave region, Pipe = pipe region
    "StaveTmin":     999., 
    "PipeTmax":    -999., 
    "PipeTmin":     999., 
  }

  def __init__ (self, roo_name, cfg_name = "config_frame", fig_outdir = "plot") :
    #
    # better run root on batch mode
    #
    ROOT.gROOT.SetBatch()

    #
    # create the output folders if not exist yet
    #
    self._fig_outdir = fig_outdir
    if not os.path.isdir( fig_outdir ):
      os.mkdir( fig_outdir )
    self._hist_outdir = self._fig_outdir + "/hist"
    if not os.path.isdir( self._hist_outdir ):
      os.mkdir( self._hist_outdir )
    self._fit_outdir = self._fig_outdir + "/fit"
    if not os.path.isdir( self._fit_outdir ):
      os.mkdir( self._fit_outdir )
  
    if not os.path.isfile( cfg_name ):
      print ("ERROR:<FRAMEANALYSIS::__INIT__> config file " + cfg_name + " not found.")
      raise Exception(" Config file error! ")

    os.system('cp '+cfg_name+' '+fig_outdir+'/'+cfg_name)
    frame_name = roo_name.split('/')[-1]
    print(frame_name) 
    os.system('cp '+roo_name+' '+fig_outdir+'/'+frame_name)

    _f_cfg = open( cfg_name, 'r')
    for line in _f_cfg:
      #
      # skip empty line or started with '#'
      #
      if ( len(line) <= 0 ) or line.startswith( '#' ) :
        continue

      item_val = line.split()
      if ( len(item_val) != 2 ): 
        print ("ERROR:<FRAMEANALYSIS::__INIT__> config items expected to be: Item Value. Not the correct style: " + line + ".")
        raise Exception(" Config file style error! ")
      item = item_val[0] 
      value = float( item_val[1] )
      print ("TEST, "+item+", value " + str(value) )
      self._parameters[ item ] = value
    print ("INFO:<FRAMEANALYSIS::__INIT__> the list of the parameters below: ")

    if (self._parameters[ "CMperPixel" ] <= 0) or (self._parameters[ "StavePixelX1" ] <= 0) or (self._parameters[ "PipePixelX1" ] <= 0):
      raise Exception(" Config parameter error! ")

    if (self._parameters[ "StavePixelX0" ] > self._parameters[ "PipePixelX0" ] > 0) or (self._parameters[ "StavePixelX1" ] < self._parameters[ "PipePixelX1" ] > 0):
      raise Exception(" Stave and pipe pixel not matched error! ")

    self._nxpixel_stave = int(self._parameters[ "StavePixelX1" ]) - int(self._parameters[ "StavePixelX0" ]) + 1
    self._nypixel_stave = int(self._parameters[ "StavePixelY1" ]) - int(self._parameters[ "StavePixelY0" ]) + 1
    self._nxpixel_pipe  = int(self._parameters[ "PipePixelX1" ] ) - int(self._parameters[ "PipePixelX0" ] ) + 1
    self._nypixel_pipe  = int(self._parameters[ "PipePixelY1" ] ) - int(self._parameters[ "PipePixelY0" ] ) + 1

    #
    # conversion from pixel to CM, pipe should use the same diameter as in Stave's coordinates
    #
    self._X0_pipe = self._parameters[ "CMperPixel" ] * ( self._parameters[ "PipePixelX0" ] - self._parameters[ "StavePixelX0" ] )
    self._Y0_pipe = self._parameters[ "CMperPixel" ] * ( self._parameters[ "PipePixelY0" ] - self._parameters[ "StavePixelY0" ] )
    self._X1_pipe = self._parameters[ "CMperPixel" ] * self._nxpixel_pipe + self._X0_pipe
    self._Y1_pipe = self._parameters[ "CMperPixel" ] * self._nypixel_pipe + self._Y0_pipe

    #
    # loop over the parameter keys
    #
    for par in self._parameters:
      val = self._parameters[ par ]
      print ("INFO:<FRAMEANALYSIS::__INIT__> " + par + " = " + str( val ) )

    #
    # read the input root file
    #
    if not os.path.isfile( roo_name ):
      print ("ERROR:<FRAMEANALYSIS::__INIT__> root input file " + roo_name + " not found.")
      raise Exception(" Root input file error! ")
    _f_roo = ROOT.TFile( roo_name, "read" );

    _btree = _f_roo.Get("btree");
    _nxpixel  = numpy.zeros(1, dtype=int)
    _nypixel  = numpy.zeros(1, dtype=int)
    _btree.SetBranchAddress( "nxpixel", _nxpixel )
    _btree.SetBranchAddress( "nypixel", _nypixel )
    _btree.GetEntry( 0 )
    if (_nxpixel[0] <=0) or (_nypixel[0] <=0):
      print ("ERROR:<FRAMEANALYSIS::__INIT__> number of pixels in X and/or Y not obtained.")
      raise Exception(" Number of pixels not set!")

    if ( self._parameters[ "StavePixelX1" ] >= _nxpixel[0] ) or ( self._parameters[ "StavePixelY1" ] >= _nypixel[0] ):
      raise Exception(" Stave pixel index overflowed error! ")
 
    self._nxpixel_raw = _nxpixel[0]
    self._nypixel_raw = _nypixel[0]

    #
    # Temperature in 2D for the whole raw figure T[y][x] initialized with -999 C.
    #
    self.stave_temperature_2d = [[ -999. for x in range( _nxpixel[0] )] for y in range( _nypixel[0] )]

    _atree = _f_roo.Get("atree");
    _xpos  = numpy.zeros(1, dtype=int)
    _ypos  = numpy.zeros(1, dtype=int)
    _temperature = numpy.zeros(1, dtype=float) 
    _atree.SetBranchAddress( "xpos", _xpos )
    _atree.SetBranchAddress( "ypos", _ypos )
    _atree.SetBranchAddress( "temperature", _temperature )
    _n_entries = _atree.GetEntries()
    for ientry in range( _n_entries ):
      _atree.GetEntry( ientry )
      if ( self._parameters[ "StaveSideL" ] ):
      # 
      # for L side, use a mirror image for Y axis to present the stave as J side
      #
        _ypos_mr = _nypixel[0] - _ypos[0] - 1
        self.stave_temperature_2d[ _ypos_mr ][ _xpos[0] ] = _temperature[0]
      else:
        self.stave_temperature_2d[ _ypos[0] ][ _xpos[0] ] = _temperature[0]

    _f_roo.Close()

  def draw_frames(self):
    """
    @brief: draw 2D Temperature frames 
      raw frame: including all pixels
      stave frame: only the stave area
      pipe frame: only the pipe area
    """

    _Xcm_frame = self._parameters[ "CMperPixel" ] * self._nxpixel_raw
    _Ycm_frame = self._parameters[ "CMperPixel" ] * self._nypixel_raw
    _h2_frame = ROOT.TH2F( "frame", ";X in cm; Y in cm; Temperature (#circC)", self._nxpixel_raw, 0., _Xcm_frame, self._nypixel_raw, 0., _Ycm_frame )

    _Xcm_stave = self._parameters[ "CMperPixel" ] * self._nxpixel_stave
    _Ycm_stave = self._parameters[ "CMperPixel" ] * self._nypixel_stave
    _h2_stave = ROOT.TH2F( "stave", ";X in cm; Y in cm; Temperature (#circC)", self._nxpixel_stave, 0., _Xcm_stave, self._nypixel_stave, 0., _Ycm_stave )

    _h2_pipe = ROOT.TH2F( "pipe", ";X in cm; Y in cm; Temperature (#circC)", self._nxpixel_pipe, self._X0_pipe, self._X1_pipe, self._nypixel_pipe, self._Y0_pipe, self._Y1_pipe )

    if ( self._parameters[ "FrameTmax" ] > -998. ):
      _h2_frame.SetMaximum( self._parameters[ "FrameTmax" ] )
    if ( self._parameters[ "FrameTmin" ] <  998. ):
      _h2_frame.SetMinimum( self._parameters[ "FrameTmin" ] )

    if ( self._parameters[ "StaveTmax" ] > -998. ):
      _h2_stave.SetMaximum( self._parameters[ "StaveTmax" ] )
    if ( self._parameters[ "StaveTmin" ] <  998. ):
      _h2_stave.SetMinimum( self._parameters[ "StaveTmin" ] )

    if ( self._parameters[ "PipeTmax" ] > -998. ):
      _h2_pipe.SetMaximum( self._parameters[ "PipeTmax" ] )
    if ( self._parameters[ "PipeTmin" ] <  998. ):
      _h2_pipe.SetMinimum( self._parameters[ "PipeTmin" ] )


    for ix in range ( self._nxpixel_raw ) :
      for iy in range ( self._nypixel_raw ) :
        T = self.stave_temperature_2d[ iy ][ ix ]
        _h2_frame.SetBinContent( ix + 1, iy + 1, T);
        if ( ix >= self._parameters[ "StavePixelX0" ] ) and ( iy >= self._parameters[ "StavePixelY0" ] ) and \
           ( ix <= self._parameters[ "StavePixelX1" ] ) and ( iy <= self._parameters[ "StavePixelY1" ] ):
          ix_stave = int( ix - self._parameters[ "StavePixelX0" ] + 1 )
          iy_stave = int( iy - self._parameters[ "StavePixelY0" ] + 1 )
          _h2_stave.SetBinContent( ix_stave, iy_stave, T)
        if ( ix >= self._parameters[ "PipePixelX0" ] ) and ( iy >= self._parameters[ "PipePixelY0" ] ) and \
           ( ix <= self._parameters[ "PipePixelX1" ] ) and ( iy <= self._parameters[ "PipePixelY1" ] ):
          ix_pipe = int( ix - self._parameters[ "PipePixelX0" ] + 1 )
          iy_pipe = int( iy - self._parameters[ "PipePixelY0" ] + 1 )
          _h2_pipe.SetBinContent( ix_pipe, iy_pipe, T)


    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas( 'c1', '', 800, 600 )
    # margin: Float_t left, Float_t right, Float_t bottom, Float_t top
    c1.SetMargin( 0.08, 0.15, 0.1, 0.05)
    _h2_frame.Draw("colz")
    c1.Print( self._fig_outdir + "/frame.png" )
    c1.Print( self._fig_outdir + "/frame.pdf" )

    c2 = ROOT.TCanvas( 'c2', '', 700, 400 )
    c2.SetMargin( 0.0675, 0.12, 0.12, 0.05)
    _h2_stave.GetXaxis().SetTitleSize( 1.3 * _h2_stave.GetXaxis().GetTitleSize() );
    _h2_stave.GetYaxis().SetTitleSize( 1.3 * _h2_stave.GetYaxis().GetTitleSize() );
    _h2_stave.GetYaxis().SetTitleOffset( 0.65 * _h2_stave.GetYaxis().GetTitleOffset() );
    _h2_stave.Draw("colz")
    c2.Print( self._fig_outdir + "/stave.png" )
    c2.Print( self._fig_outdir + "/stave.pdf" )

    #c3 = ROOT.TCanvas( 'c3', '', 700, 400 )
    #c3.SetMargin( 0.0675, 0.12, 0.12, 0.05)
    _h2_pipe.GetXaxis().SetTitleSize( 1.3 * _h2_pipe.GetXaxis().GetTitleSize() );
    _h2_pipe.GetYaxis().SetTitleSize( 1.3 * _h2_pipe.GetYaxis().GetTitleSize() );
    _h2_pipe.GetYaxis().SetTitleOffset( 0.65 * _h2_pipe.GetYaxis().GetTitleOffset() );
    _h2_pipe.Draw("colz")
    c2.Print( self._fig_outdir + "/pipe.png" )
    c2.Print( self._fig_outdir + "/pipe.pdf" )

  def fit_hist(self, h1):
    """
      read a 1D histogram, fit to gaussian function, return (mean, mean position, width) 
      expecting h1 in range [0., n] with n the number of bins
    """
   
    ROOT.gStyle.SetOptStat(0)
    cf1 = ROOT.TCanvas( 'cf1', '', 600, 600 )
    # margin: Float_t left, Float_t right, Float_t bottom, Float_t top
    cf1.SetMargin( 0.1, 0.05, 0.12, 0.05)
    h1.GetXaxis().SetTitle( "Y axis in pixels" )
    h1.GetYaxis().SetTitle( "Temperature (#circC)")
    h1.Draw()

    nbins = h1.GetNbinsX()
    f1 = ROOT.TF1( "FitGaus", "gaus",0.0, float(nbins) ) 
    if ( self._parameters[ "LiquidTLow" ] > 0 ):
      f1.SetParLimits(0, -100., 10.)
    else:
      f1.SetParLimits(0, 0.00001, 100.);

    #constraints on sigma
    f1.SetParameter(2, nbins/3. ); 
    f1.SetParLimits(2, nbins/5., nbins * 3 / 4.); 

    #constraints on peak position
    f1.SetParameter(1, nbins/2.);
    f1.SetParLimits(1, nbins/4., nbins * 3 / 4.);
    h1.Fit(f1,"Q") # Fit(f1,"0") meaning no draw on the canvas

    # 3 parameters: temperature, mean, sigma
    _pars  = numpy.zeros(3, dtype=float)
    f1.GetParameters( _pars );
 
    _temp = _pars[0]
    _mean = _pars[1]
    _width= _pars[2]
    _chi2 = f1.GetChisquare()
    _ndf  = f1.GetNDF()

    Tl = ROOT.TLatex()
    Tl.SetTextSize(20)
    Tl.SetTextFont(43)
    Tl.SetNDC()
    ss = "#color[4]{Left: #chi^{2}/NDF = %.1f / %d}" % (_chi2, _ndf)
    Tl.DrawLatex(0.18, 0.88, ss)
    Tl.DrawLatex(0.18, 0.84, ("#color[4]{mean: %5.3f}" % _mean))
    Tl.DrawLatex(0.18, 0.80, ("#color[4]{width: %5.3f}" % _width))
    Tl.DrawLatex(0.18, 0.76, ("#color[4]{T: %3.1f #circ C}" % _temp))

    cf1.Update()

    _name = self._fit_outdir + "/h_" + str( h1.GetTitle() ) + ".png"
    cf1.Print( _name )
 
    return (_temp, _mean, _width, _chi2, _ndf)

  def find_pipes(self):
    """
      Find the temperature profile of the cooling pipe. At each point along the cooling pipe, find the minimum ( or maximum depending on the operating
        temperature ) along Y axis.

        Cooling pipe looks like:
         Y -----------
                      
           -----------
                       X ->
        where, only the top and bottom pipes are kept. The short pipe turn along the Y axis at right is NOT analyzed so far.
        Obtain the cooling pipe temperature as a function of X axis for top and bottom lines. Use the top 40% of the pixels to get the top cooling pipe,
        and the bottom 40% of the pixels for the bottom pipe curve.
    """

    _roo_out = ROOT.TFile("plot/result.root", "recreate")
    h1s = [ ROOT.TH1F("top_pipe_temperature",  ";X in cm; Temperature (#circC)",  self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("top_pipe_mean",         ";X in cm; Pipe position Y in cm", self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("top_pipe_width",        ";X in cm; Pipe width Y in cm",    self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("top_pipe_chi2",         ";X in cm; Pipe fit chi2",         self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("top_pipe_ndf",          ";X in cm; Pipe fit ndf",          self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("bot_pipe_temperature",  ";X in cm; Temperature (#circC)",  self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("bot_pipe_mean",         ";X in cm; Pipe position Y in cm", self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("bot_pipe_width",        ";X in cm; Pipe width Y in cm",    self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("bot_pipe_chi2",         ";X in cm; Pipe fit chi2",         self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
            ROOT.TH1F("bot_pipe_ndf",          ";X in cm; Pipe fit ndf",          self._nxpixel_pipe, self._X0_pipe, self._X1_pipe),
          ]
          

    for ix in range ( self._nxpixel_pipe ) :
      ix_raw = int ( ix + self._parameters[ "PipePixelX0" ] )
      ny_1pipe = int ( 0.4 * self._nypixel_pipe )
      #
      # pipe using the top and bottom 40% of the pixels
      #   top: high Y
      #   bottom: low Y
      #
      ht1 = ROOT.TH1F( "t"+str(ix), "t"+str(ix), ny_1pipe, 0., float( ny_1pipe ) ) 
      hb1 = ROOT.TH1F( "b"+str(ix), "b"+str(ix), ny_1pipe, 0., float( ny_1pipe ) ) 
      for iy in range ( self._nypixel_pipe ) :
        
        iy_raw = int ( iy + self._parameters[ "PipePixelY0" ] )

        if ( iy < ny_1pipe ):
          #
          # low Y pixel number for bottom curve
          #
          hb1.SetBinContent( iy + 1, self.stave_temperature_2d[ iy_raw ][ ix_raw ])
          hb1.SetBinError( iy + 1, 0.02 * self.stave_temperature_2d[ iy_raw ][ ix_raw ]) # set 2% error
        elif ( iy >= int(self._nypixel_pipe - ny_1pipe) ):
          #
          # high Y pixel number for top curve
          #
          iy_reset = int(iy - (self._nypixel_pipe - ny_1pipe) + 1 )
          ht1.SetBinContent( iy_reset, self.stave_temperature_2d[ iy_raw ][ ix_raw ])
          ht1.SetBinError( iy_reset, 0.02 * self.stave_temperature_2d[ iy_raw ][ ix_raw ]) # set 2% error

      # returned tuple: (temp, mean, width, chi2, ndf)
      _t_data = self.fit_hist( ht1 )
      _b_data = self.fit_hist( hb1 )
      for idx,val in enumerate(_t_data):
        h1s[ idx ].SetBinContent( ix+1, val )
      for idx,val in enumerate(_b_data):
        jdx = int( idx + 5 )
        h1s[ jdx ].SetBinContent( ix+1, val )

    c0 = ROOT.TCanvas( 'c0', '', 2000, 600 )
    # margin: Float_t left, Float_t right, Float_t bottom, Float_t top
    c0.SetMargin( 0.065, 0.03, 0.15, 0.05)
 
    for h1 in h1s:
      h1.SetLineWidth(3)
      h1.GetXaxis().SetTitleSize( 1.3 * h1.GetXaxis().GetTitleSize() );
      h1.GetYaxis().SetTitleSize( 1.3 * h1.GetYaxis().GetTitleSize() );
      h1.GetYaxis().SetTitleOffset( 0.65 * h1.GetYaxis().GetTitleOffset() );
      h1.Draw()
      c0.Print( self._hist_outdir + "/" + h1.GetName() + ".png" )
      c0.Print( self._hist_outdir + "/" + h1.GetName() + ".pdf" )

      h1.Write()

    _roo_out.Close()

def print_usage( s_function):
  print ("Usage: " + s_function + " INPUT_ROOT_FILE [CONFIG=config_frame] [OUT_DIR=plot]")

def main():
  if sys.version_info[0] >= 3:
    print ("ERROR:<FRAMEANALYSIS::MAIN> PyROOT only works with Python 2.x. Code Tested with 2.7.10. Current version " + str(sys.version_info[0]) + ".x")
    raise Exception(" Python Version too high. Use 2.x. ")

  nargv = len(sys.argv)

  argv0 = str(sys.argv[0])

  bolFindConfig = True
  if (nargv <= 1): 
    print ("ERROR:<FRAMEANALYSIS> Please provide: input root file. Missing! Return.")
    print_usage( argv0 )
    return
  elif ( argv0 == '-h' ) or ( argv0 == '--help' ):
    print_usage( argv0 )
    return
  elif ( sys.argv[1] == '-mc') or ( sys.argv[1] == '--manualconfig' ):
    print ("Usage: Manual configuration overide. Will use config settings from old config_frame")
    bolFindConfig = False
    str_inroo = sys.argv[2];
  else:
    str_inroo = sys.argv[1];

  str_cfg = "config_frame"
  if (nargv >= 3) and bolFindConfig == True:
    str_cfg = sys.argv[2];
  elif (nargv >= 4) and bolFindConfig == False:
    str_cfg = sys.argv[3];

  str_outdir = "plot"
  if (nargv >= 4) and bolFindConfig == True:
    str_outdir = sys.argv[3];
  elif (nargv >= 5) and bolFindConfig == False:
    str_outdir = sys.argv[4];

  if bolFindConfig == True:
    print ("Usage: Finding frame configuration...")
    cf.FindPoints( str_inroo, str_cfg, str_outdir)

  ist_frmana = FrameAnalysis( str_inroo, str_cfg, str_outdir )
  ist_frmana.draw_frames()
  ist_frmana.find_pipes()
  print (' Make plots. Done!')

if __name__ == "__main__":
  main()
