#!/usr/bin/python
"""
@run:
  ./seqToProfile.py [files]

  [files]: This can be any number of differently named .seq files

@brief:
  The code will take each seq file and convert it to a root file using
  read_sequence.sh. It will then run frameanal.py on the frame_average.root
  file created by read_sequence.sh. frameanal.py initially will put
  everything into a folder called plots. This will then be relabeled by the
  same name as the sequence file.
"""

import sys
import os
import ROOT
import numpy as np
import time

def main():
  """
    The main loop
  """
  #load the files
  nargv = len(sys.argv)
  inputfiles = []
  if (nargv <=1):
    print("ERROR: Please provied a set of files.")
  else:
    for i in range(1,nargv):
      inputfiles = np.append(inputfiles,sys.argv[i])
  while True:
    Vin = raw_input("\nIs the stave core a 13 or 14 module core? (13/14)")
    if "13" in Vin:
      bol14Mod = False
      break
    elif "14" in Vin:
      bol14Mod = True
      break
    else:
      print("\n Please provide either a 13 or 14")

  starttime = time.strftime('%a, %d %b %Y %H:%M:%S',time.localtime())
  print('Start Time: '+str(starttime))
  print(str(len(inputfiles))+' files will be converted...')
  print(inputfiles)

  #Do the thing!
  nfiles = len(inputfiles)
  for i in range(nfiles):
    ctime = time.strftime('%a, %d %b %Y %H:%M:%S',time.localtime())
    print('BEGINNING LOOP {0}/{1}'.format(i+1,nfiles))
    print('  Start Time   : '+str(starttime))  
    print('  Current Time : '+str(ctime))
    #Read the sequence
    os.system('bash read_sequence.sh '+inputfiles[i]+' -e 0.92')
    #Create the plots
    if bol14Mod == False:
      os.system('./frameanal.py roo/frame_average.root')
    else:
      os.system('./frameanal.py roo/frame_average.root -14M')
    #Rename the outdir
    outfilename = inputfiles[i].split('.')[0]
    try:
      os.system('mv plot plot-'+outfilename)
    except:
      os.system('mv roo/frame_average.root '+outfilename)
  endtime = time.strftime('%a, %d %b %Y %H:%M:%S',time.localtime())

  print('Start Time: '+str(starttime))  
  print('End Time  : '+str(endtime))


if __name__ == '__main__':
  main()
