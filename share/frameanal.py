#!/usr/bin/python

"""
@run

@brief:
  
@functions (class TextToRoot):

@notes:
  This code assumes using python 2.7.10
  It doesn't work for python 3.x, because ROOT support is not guaranteed.

@reference:
  TTree in python: 
    https://www-zeuthen.desy.de/~middell/public/pyroot/pyroot.html
    https://root.cern.ch/how/how-write-ttree-python

@email: jie.yu@cern.ch
"""

import sys
import os
import numpy
import math

# running root in batch mode
# see: https://root-forum.cern.ch/t/pyroot-in-batch-mode/3105/2
# sys.argv.append( '-b-' )
import ROOT

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
    ROOT.gROOT.SetBatch()

    self._fig_outdir = fig_outdir
    if not os.path.isdir( fig_outdir ):
      os.mkdir( fig_outdir )
    self._hist_outdir = self._fig_outdir + "/hist"
    if not os.path.isdir( self._hist_outdir ):
      os.mkdir( self._hist_outdir )
 
    if not os.path.isfile( cfg_name ):
      print ("ERROR:<FRAMEANALYSIS::__INIT__> config file " + cfg_name + " not found.")
      raise Exception(" Config file error! ")
    _f_cfg = open( cfg_name, 'r')
    for line in _f_cfg:
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

    #
    # Pipe should use the same diameter as in Stave's coordinates
    #
    _X0_pipe = self._parameters[ "CMperPixel" ] * ( self._parameters[ "PipePixelX0" ] - self._parameters[ "StavePixelX0" ] )
    _Y0_pipe = self._parameters[ "CMperPixel" ] * ( self._parameters[ "PipePixelY0" ] - self._parameters[ "StavePixelY0" ] )
    _X1_pipe = self._parameters[ "CMperPixel" ] * self._nxpixel_pipe + _Y0_pipe
    _Y1_pipe = self._parameters[ "CMperPixel" ] * self._nypixel_pipe + _Y0_pipe

    _h2_pipe = ROOT.TH2F( "pipe", ";X in cm; Y in cm; Temperature (#circC)", self._nxpixel_pipe, _X0_pipe, _X1_pipe, self._nypixel_pipe, _Y0_pipe, _Y1_pipe )

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
    f1 = ROOT.TF1( "FitGaus", "gaus", 0., float(nbins) ) 
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
    h1.Fit(f1) # Fit(f1,"0") meaning no draw on the canvas

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

    _name = self._hist_outdir + "/h_" + str( h1.GetTitle() ) + ".png"
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

    for ix in range ( self._nxpixel_pipe ) :
      ix_raw = int ( ix + self._parameters[ "PipePixelX0" ] )
      ny_1pipe = int ( 0.4 * self._nypixel_pipe )
      #
      # pipe using the top and bottom 40% of the pixels
      #
      ht1 = ROOT.TH1F( "t"+str(ix), "t"+str(ix), ny_1pipe, 0., float( ny_1pipe ) ) 
      hb1 = ROOT.TH1F( "b"+str(ix), "b"+str(ix), ny_1pipe, 0., float( ny_1pipe ) ) 
      for iy in range ( self._nypixel_pipe ) :
        
        iy_raw = int ( iy + self._parameters[ "PipePixelY0" ] )

        if ( iy < ny_1pipe ):
          ht1.SetBinContent( iy + 1, self.stave_temperature_2d[ iy_raw ][ ix_raw ])
          ht1.SetBinError( iy + 1, 0.02 * self.stave_temperature_2d[ iy_raw ][ ix_raw ]) # set 2% error
        elif ( iy >= int(self._nypixel_pipe - ny_1pipe) ):
          iy_reset = int(iy - (self._nypixel_pipe - ny_1pipe) + 1 )
          hb1.SetBinContent( iy_reset, self.stave_temperature_2d[ iy_raw ][ ix_raw ])
          hb1.SetBinError( iy_reset, 0.02 * self.stave_temperature_2d[ iy_raw ][ ix_raw ]) # set 2% error
      _t_temp, _t_mean, _t_width, _t_chi2, _t_ndf = self.fit_hist( ht1 )
      _b_temp, _b_mean, _b_width, _b_chi2, _b_ndf = self.fit_hist( hb1 )

def print_usage( s_function):
  print ("Usage: " + s_function + " INPUT_ROOT_FILE [CONFIG=config_frame] [OUT_DIR=plot]")

def main():
  if sys.version_info[0] >= 3:
    print ("ERROR:<FRAMEANALYSIS::MAIN> PyROOT only works with Python 2.x. Code Tested with 2.7.10. Current version " + str(sys.version_info[0]) + ".x")
    raise Exception(" Python Version too high. Use 2.x. ")

  nargv = len(sys.argv)
  if (nargv <= 1): 
    print ("ERROR:<FRAMEANALYSIS> Please provide: input root file. Missing! Return.")
    print_usage( str(sys.argv[0]) )
    return
  else:
    str_inroo = sys.argv[1];

  str_cfg = "config_frame"
  if (nargv >= 3):
    str_cfg = sys.argv[2];

  str_outdir = "plot"
  if (nargv >= 4):
    str_outdir = sys.argv[3];

  ist_frmana = FrameAnalysis( str_inroo, str_cfg, str_outdir )
  ist_frmana.draw_frames()
  ist_frmana.find_pipes()
  print (' Make plots. Done!')

if __name__ == "__main__":
  main()
