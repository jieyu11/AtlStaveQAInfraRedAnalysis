#!/usr/bin/env python
# @purpose:
#   Read .seq files from Flir IR camera and write each frame to temporary binary file.
#
# @usage:
#   seqToBin.py _FILE_NAME_.seq
#
# @note:
#   When first using this code for a new camera, it might need find the bits separating
#   each frame, which is possibly IR camera specific. Please run:
#     hexdump -n16 -C _FILE_NAME_.seq 
#
#   @@Example
#     >$ hexdump -n16 -C Rec-000667_test.seq 
#     00000000  46 46 46 00 52 65 73 65  61 72 63 68 49 52 00 00  |FFF.ResearchIR..|
#     00000010
#     So, for this camera, the separation patten is:
#     \x46\x46\x46\x00\x52\x65\x73\x65\x61\x72\x63\x68\x49\x52
#     which == FFFResearchIR
#

import sys
import os

#pat=b'\x46\x46\x46\x00\x52\x65\x73\x65\x61\x72\x63\x68\x49\x52';
pat='\x46\x46\x46\x00\x52\x65\x73\x65\x61\x72\x63\x68\x49\x52'
#pat = 'FFF.ResearchIR'

def split_by_marker(f, marker = pat, block_size = 10240):
  current = ''
  bolStartPos = True
  while True:
    block = f.read(block_size)
    if not block: # end-of-file
      yield marker+current
      return
    current += block
    while True:
      markerpos = current.find(marker) 
      if bolStartPos ==True:
        current = current[markerpos +len(marker):]
        bolStartPos = False
        continue
      elif markerpos <0:
        break
      else:
        yield marker+current[:markerpos]
        current = current[markerpos+ len(marker):]

def main():
  print "This is the name of the script: ", sys.argv[0]
  print "Number of arguments: ", len(sys.argv)
  print "The arguments are: " , str(sys.argv)

  idx=0
  inputname=sys.argv[1]
  outdir = sys.argv[2]
  if not os.path.isdir(outdir):
    os.mkdir(outdir)
  for line in split_by_marker(open(inputname, 'rb')):
    outname="fout/frame_{0}.fff".format( idx )
    output_file = open( outname ,"wb")
    output_file.write( line )
    output_file.close()
    print outname
    idx=idx+1
    if idx % 100000 == 0:
      print 'running index : {} '.format( idx )
      break;


if __name__=='__main__':
  main()
