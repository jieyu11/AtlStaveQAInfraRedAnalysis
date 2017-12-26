* Purpose:
  
  - Convert raw ADC counts recorded by the Flir InfraRed Camera into temperature (C).
  - Analyze the stave temperature profile and extract the temperature profile of the cooling pipe in the stave.

* Run:
  - Check installation below first!

  - ./read_sequence.sh FlirIR_rec_example.seq
    + convert .seq to .root 
    + each image is stored in one root file roo/frame_idx.root, idx = [0, ...]
    + an average image is stored in file roo/frame_average.root

  - ./frameanal.py roo/frame_average.root
    + plots of frames in plot/.
    + plots of cooling pipe temperature profile in plot/hist/.
    + plots of fit temperature to get min/max around cooling pipe in plot/fit/.


* Installation:

  - There are a few softwares needed before running this package. A list below.

  - ExifTool (read .seq file produced by Flir IR camera)
    + download: http://www.sno.phy.queensu.ca/~phil/exiftool/
    + information: http://www.sno.phy.queensu.ca/~phil/exiftool/TagNames/FLIR.html (e.g.: .seq to .png)
    + command: "exiftool -Emissivity Test.seq" gives "Emissivity : 0.95"
    + command: "exiftool -FLIR:all Test.seq" gives all properties related to the .seq
    + command: "exiftool -RawThermalImage -b  Test.seq > firstImage.dat" extracting the first image from Test.seq

  - MacPorts (needed for command 'convert', see share/binarytotext.sh)  
    + download: https://www.macports.org/install.php

  - ImageMagick (needed for command 'convert') 
    + download: http://cactuslab.com/imagemagick/

