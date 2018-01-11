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
    + for help, use: "./frameanal.py -h (OR --help)"
    + modify 'config_frame' correspondingly to your setup
      1. Pixel number for the Stave region
      2. Pixel number for the Pipe region
      3. centimeter per pixel
      4. is L side
      5. is temperature lower than room temperature
      6. maximum and mimimum temperature for plots of Frame, Stave and Pipe
    + out plots of frames in plot/.
    + out plots of cooling pipe temperature profile in plot/hist/.
    + out plots of fit temperature to get min/max around cooling pipe in plot/fit/.
    + out cooling pipe result also stored in result.root


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

  - ROOT (from CERN to make plots and store data)
    + download: https://root.cern.ch/downloading-root
