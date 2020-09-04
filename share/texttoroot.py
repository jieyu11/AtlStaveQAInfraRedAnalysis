#!/usr/bin/env python 

"""
@run
  ./share/texttoroot.py OUT_DIR NUM_INPUT_FILES [CONFIG=config] [IN_DIR=tout] [IN_NAME=frame] [IN_EXT=pgm]

  parameters:
    OUT_DIR: output directory, necessary
    NUM_INPUT_FILES: number of input files, necessary
    CONFIG: IR camera configuration file, optional, default: config
    IN_DIR: input file directory, optional, default: tout
    IN_NAME: input file name prefix, optional, default: frame
    IN_EXT: input file extension, optional, default: pgm

@brief:
  This code converts ADC counts recorded by IR camera into temperature values
  in degree C.

  The text files in the following format with an example of the first 5 lines:
    1 Time 2017:08:28 13:22:56.744
    2 P2
    3 640 480
    4 65535
    5 13432 13431 13436 13461 13433 
    6 ... ...
  The first line lists the time information when the frame is recorded.
  The third line indicates the number of pixels in X direction and Y direction.
  Starting from the fifth line, the ADC counts from [x,y] of [0,0], [1,0] ...
 
  The formula used to convert ADC counts into temperature (C) is defined in
  function: counts_to_temperature( counts ) return temperature.
  Search for Planck's law to get more information.
    https://en.wikipedia.org/wiki/Planck%27s_law

  Two trees created in the root file:
  TTree: atree
  TBranch: xpos: X pixel index [0, NXPIXEL )
  TBranch: ypos: Y pixel index [0, NYPIXEL )
  TBranch: temperature at (xpos, ypos)

  - Note that the Flir IR camera starting Y axis pixels from top to bottom
    Example of [x,y]:
    [0,0], [1,0], ..., [639,0]
    [0,1], [1,1], ..., [639,1]
    ...
    But, a graph in root starting the Y axis from bottom. So, if one plots the 
    temperature as a function of [x,y] position in root, the figure is upside down.
    To avoid this, the Y axis number stored in root is reverted in this code.

  TTree: btree
  Frame time information and number of pixels in X and Y directions.

  A root file with average temperature frame is calculated and kept.
  
@functions (class TextToRoot):
  counts_to_temperature( counts ) return temperature
  - convert ADC counts into temperature pixel by pixel

  convert(outdir, n_inputs, indir = "tout", inname = "frame", inext = "pmg")
  - set the output directory and the number of input text files
  - convert each text file (representing one frame) into a root file.
  - default input directory and input file name and extensions, the input files 
    are assumed:
      tout/frame_0.pmg
      tout/frame_1.pmg
      ... ...
  - the output files are stored, with one with average temperature:
      $outdir/frame_0.root
      $outdir/frame_1.root
      ... ...
      $outdir/frame_average.root


@notes:
  This code assumes using python 2.7.10
  It doesn't work for python 3.x, because ROOT support is not guaranteed.

@reference:
  TTree in python: 
    https://www-zeuthen.desy.de/~middell/public/pyroot/pyroot.html
    https://root.cern.ch/how/how-write-ttree-python

@email: jie.yu@cern.ch
"""

import ROOT
import sys
import os
import numpy
import math
import resource
import gc

class TextToRoot:
  """

  """
  _parameters = { "R1": 0., "R2": 0., "B": 0., "O": 0., "F": 0., "Emissivity": 0.79, "ReflTemp": 20., "AtomTemp": 20., "Transmissivity": 0.987}
  #Emissivity of the stave is 0.9, pipefoam is 0.79
  def __init__ (self, cfg_name = "config") :
    self._status = 0
    if not os.path.isfile( cfg_name ):
      print ("ERROR:<TEXTTOROOT::__INIT__> config file " + cfg_name + " not found. Status = 1. ")
      _status = 1;
    _f_cfg = open( cfg_name, 'r')
    for line in _f_cfg:
      item_val = line.split()
      if ( len(item_val) != 2 ): 
        print ("ERROR:<TEXTTOROOT::__INIT__> config items expected to be: Item Value. Not the correct style: " + line + ". Status = 2.")
        _status = 2
      item = item_val[0] 
      value = float( item_val[1] )
      self._parameters[ item ] = value
    print ("INFO:<TEXTTOROOT::__INIT__> the list of the parameters below: ")

    #
    # loop over the parameter keys
    #
    for par in self._parameters:
      val = self._parameters[ par ]
      if ( math.fabs(val) < 1.e-5 ):
        print ("ERROR:<TEXTTOROOT::__INIT__> " + par + " = " + str( val ) + ". Please set! Status = 3.")
        _status = 3
      else:
        print ("INFO:<TEXTTOROOT::__INIT__> " + par + " = " + str( val ) )
    

  def get_status(self) :
    return self._status

  def counts_to_temperature(self, raw_counts ):
    """
    @brief: converting DC counts received by each pixel of the IR camera,
       which represents energy (heat), to temperature in degree C.
    """
    if (self._status > 0):
      return -999.
  
    RawAtom = self._parameters["R1"] / (self._parameters["R2"] * ( math.exp( self._parameters["B"] / (self._parameters["AtomTemp"] + 273.15) ) - self._parameters["F"] ) ) - self._parameters["O"];
    RawRefl = self._parameters["R1"] / (self._parameters["R2"] * ( math.exp( self._parameters["B"] / (self._parameters["ReflTemp"] + 273.15) ) - self._parameters["F"] ) ) - self._parameters["O"];
    RawObj_numerator   = ( raw_counts - self._parameters["Transmissivity"] * (1 - self._parameters["Emissivity"]) * RawAtom - (1 - self._parameters["Transmissivity"]) * RawRefl ) 
    RawObj_denominator = ( self._parameters["Emissivity"] * self._parameters["Transmissivity"]);
    RawObj = RawObj_numerator / RawObj_denominator
    TinC = self._parameters["B"] / math.log ( self._parameters["R1"] / ( self._parameters["R2"] * ( RawObj + self._parameters["O"] ) ) + self._parameters["F"] ) - 273.15;
  
    return TinC;

  def convert(self, outdir, n_inputs, indir = "tout", inname = "frame", inext = "pmg"):
    """
    """
    if ( n_inputs <= 0 ):
      print ("ERROR:<TEXTTOROOT::CONVERT> no input in folder " + indir + "! Return!")
      return
  
    #
    # average frame: 
    #   using 1st frame information of time
    #   using average temperature
    #

    #
    # creating 1D arrays to keep the data!!
    #
    avg_nxpixel = numpy.zeros(1, dtype=int)
    avg_nypixel = numpy.zeros(1, dtype=int)
    avg_year = numpy.zeros(1, dtype=int)
    avg_month = numpy.zeros(1, dtype=int)
    avg_date = numpy.zeros(1, dtype=int)
    avg_hour = numpy.zeros(1, dtype=int)
    avg_minute = numpy.zeros(1, dtype=int)
    avg_second = numpy.zeros(1, dtype=float)
  
    avg_temperature = numpy.zeros(1, dtype=float)
    avg_xpos  = numpy.zeros(1, dtype=int)
    avg_ypos  = numpy.zeros(1, dtype=int)
   
    strRooName_avg = outdir + "/" + inname + "_average.root" 
    f_roo_avg = ROOT.TFile( strRooName_avg, "recreate")
  
    avg_atree = ROOT.TTree("atree", "a tree of temperature data");
    avg_atree.Branch('temperature', avg_temperature, 'temperature/D')
    avg_atree.Branch('xpos', avg_xpos, 'xpos/I')
    avg_atree.Branch('ypos', avg_ypos, 'ypos/I')
  
    avg_btree = ROOT.TTree("btree", "a tree of camera information");
    avg_btree.Branch('nxpixel', avg_nxpixel, 'nxpixel/I')
    avg_btree.Branch('nypixel', avg_nypixel, 'nypixel/I')
    avg_btree.Branch('year', avg_year, 'year/I')
    avg_btree.Branch('month', avg_month, 'month/I')
    avg_btree.Branch('date', avg_date, 'date/I')
    avg_btree.Branch('hour', avg_hour, 'hour/I')
    avg_btree.Branch('minute', avg_minute, 'minute/I')
    avg_btree.Branch('second', avg_second, 'second/D') 

    #Things that didn't need to be in the loop
    temperature = numpy.zeros(1, dtype=float)
    xpos  = numpy.zeros(1, dtype=int)
    ypos  = numpy.zeros(1, dtype=int)
    index  = numpy.zeros(1, dtype=int)
    nxpixel = numpy.zeros(1, dtype=int)
    nypixel = numpy.zeros(1, dtype=int)
    year = numpy.zeros(1, dtype=int)
    month = numpy.zeros(1, dtype=int)
    date = numpy.zeros(1, dtype=int)
    hour = numpy.zeros(1, dtype=int)
    minute = numpy.zeros(1, dtype=int)
    second = numpy.zeros(1, dtype=float)


    for outidx in range(n_inputs):
      #print 'Memory usage BOL: %s (MB)' % str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1048576)
      fname = indir + "/" + inname + "_" + str(outidx) + "." + inext
      if not os.path.isfile( fname ):
        print ("ERROR:<TEXTTOROOT::CONVERT> file name " + fname + " is incorrect. ")
        continue;
      else:
        print ("INFO:<TEXTTOROOT::CONVERT> converting " + fname + " to root file. ")
 
      outnum = str(outidx)
      if outidx < 10:
        outnum = "0000" + outnum
      elif outidx < 100:
        outnum = "000" + outnum
      elif outidx < 1000:
        outnum = "00" + outnum
      elif outidx < 10000:
        outnum = "0" + outnum
       
      strRooName = outdir + "/" + inname + "_" + outnum + ".root" 
      f_roo = ROOT.TFile( strRooName, "recreate")
 
      atree = ROOT.TTree("atree", "a tree of temperature data");
      atree.Branch('temperature', temperature, 'temperature/D')
      atree.Branch('xpos', xpos, 'xpos/I')
      atree.Branch('ypos', ypos, 'ypos/I')
  
      btree = ROOT.TTree("btree", "a tree of camera information");
      btree.Branch('index', index, 'index/I')
      btree.Branch('nxpixel', nxpixel, 'nxpixel/I')
      btree.Branch('nypixel', nypixel, 'nypixel/I')
      btree.Branch('index', index, 'index/I')
      btree.Branch('year', year, 'year/I')
      btree.Branch('month', month, 'month/I')
      btree.Branch('date', date, 'date/I')
      btree.Branch('hour', hour, 'hour/I')
      btree.Branch('minute', minute, 'minute/I')
      btree.Branch('second', second, 'second/D')
  
      f = open( fname, 'r')
      il = -1  # ith line
      ipix = 0 # ith pixel
      nempty_line = 0 # numbers of empty lines (shouldn't be any.)
     
      for line in f:
        il = il + 1
  
        content = line.split()
        if ( len(content) <= 0 ):
          print ("WARNING:<TEXTTOROOT::CONVERT> file " + fname + " has an empty line! Continue. ")
          nempty_line = nempty_line + 1
          continue
        if ( nempty_line >= 10 ):
          print ("ERROR:<TEXTTOROOT::CONVERT> file " + fname + " finds >=10 empty lines! Break. ")
          break
  
        #
        # Text file content example:
        #
        # Time 2017:08:28 13:22:56.744
        # P2
        # 640 480
        # 65535
        # 13432 13431 13436 13461 13433 
        #
  
        if ( il == 0 ): 
          if ( content[0] != "Time" ):
            print ("ERROR:<TEXTTOROOT::CONVERT> first line does not start with Time. Check! ")
            break
          str_ymd = content[1].split(':');
          if ( len(str_ymd) < 3 ):
            print ("ERROR:<TEXTTOROOT::CONVERT> year,month,date not all found. Check! ")
            break

          #
          # since it year itself is an array, so use year[0] to get its value
          #
          year[0] = int(str_ymd[0])
          month[0] = int(str_ymd[1]) 
          date[0] = int(str_ymd[2])
  
          str_hms = content[2].split(':');
          if ( len(str_hms) < 3 ):
            print ("ERROR:<TEXTTOROOT::CONVERT> hour,minute,date not all found. Check! ")
            break
          hour[0] = int(str_hms[0])
          minute[0] = int(str_hms[1]) 
          second[0] = float(str_hms[2])
        elif ( il == 1 ) or ( il == 3 ):  
          continue
        elif ( il == 2 ):
          if ( len(content) < 2 ):
            print ("ERROR:<TEXTTOROOT::CONVERT> number of X and Y pixels not found. Check! ")
            break
          nxpixel[0] = int( content[0] )
          nypixel[0] = int( content[1] )
          if ( outidx == 0 ):
            print ("INFO:<TEXTTOROOT::CONVERT> NXPIX " +str(nxpixel[0]) + " NYPIX " + str(nypixel[0]) )
          btree.Fill()
  
          #
          # Initialize the values with 0 for the average frame when reading the first frame
          #
          if ( outidx == 0 ):
            avg_temperature_2d = [[ 0. for x in range( nxpixel[0] )] for y in range( nypixel[0] )]
  
          #
          # keep the first frame time information for the average one!
          #
          
            avg_nxpixel[0] = nxpixel[0]
            avg_nypixel[0] = nypixel[0]
            avg_year[0] = year[0]
            avg_month[0] = month[0]
            avg_date[0] = date[0]
            avg_hour[0] = hour[0]
            avg_minute[0] = minute[0]
            avg_second[0] = second[0]
        else: 
          #
          # content is a list of raw counts
          #----------------------------MEMORY LEAK BELOW
          
          for str_count in content: 

            xpos[0] = int ( ipix % nxpixel[0] )

            #
            # note the Y axis pixel index is reverted top <--> bottom
            #
            ypos[0] = nypixel[0] - int ( ipix / nxpixel[0] ) - 1
            counts = int( str_count )
            temperature[0] = self.counts_to_temperature( counts )
            if (ipix == 0) and (outidx == 0 ):
              print ( "INFO:<TEXTTOROOT::CONVERT> ix: " + str(xpos[0]) + " iy: " +str(ypos[0]) + " count: " + str(counts) + " T: " + str(temperature[0]) )

            avg_temperature_2d[ypos[0]][xpos[0]] += temperature[0] / n_inputs
  
            ipix = ipix + 1
            atree.Fill()
          #----------------------------MEMORY LEAK ABOVE
          # end of one line
          # ---------------

        # check the line number
        # ---------------
     
      # end of all lines
      # ----------------
      #print 'Memory usage B4C: %s (MB)' % str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1048576)

      f.close()
      f_roo.Write()
      #print(f_roo.GetSize())
      
      #print (str(atree.GetTotBytes()/1048576)) 
      #print 'Memory usage EOL: %s (MB)' % str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1048576)
      
    #
    # now deal with the average
    #
    avg_btree.Fill();
    for x in range( avg_nxpixel[0] ):
      for y in range( avg_nypixel[0] ):
        avg_xpos[0] = x
        avg_ypos[0] = y
        avg_temperature[0] = avg_temperature_2d[ y ][ x ]
        avg_atree.Fill()
    f_roo_avg.Write()
    f_roo_avg.Close()
    #print 'Memory usage Fin: %s (MB)' % str(float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1048576)
def print_usage( s_function):
  print ("Usage: " + s_function + " OUT_DIR NUM_INPUT_FILES [CONFIG=config] [IN_DIR=tout] [IN_NAME=frame] [IN_EXT=pgm]")

def main():

  nargv = len(sys.argv)
  if (nargv <= 2): 
    print ("ERROR:<TEXTTOROOT> Please provide: output folder and number of inputs. Missing! Return.")
    print_usage( str(sys.argv[0]) )
    return
  else:
    str_outdir = sys.argv[1];
    int_ninput = int( sys.argv[2] )#+3214;

  str_cfg = "config"
  if (nargv >= 4):
    str_cfg = sys.argv[3];

  str_indir = "tout"
  if (nargv >= 5):
    str_indir = sys.argv[4];

  str_inname = "frame"
  if (nargv >= 6): 
    str_inname = sys.argv[5];

  str_inext = "pgm"
  if (nargv >= 7):
    str_inext = sys.argv[6];

  ist_txtroo = TextToRoot( str_cfg ) 
  ist_txtroo.convert( str_outdir, int_ninput, str_indir, str_inname, str_inext)
  print (' Convert. Done!')

if __name__ == "__main__":
  main()
