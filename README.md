* Purpose:
  
  - Convert raw ADC counts recorded by the Flir InfraRed Camera into temperature (C).
  - Analyze the stave temperature profile and extract the temperature profile of the cooling pipe in the stave.

* Run:
  - Check installation below first!

  - ./read_sequence.sh FlirIR_rec_example.seq  [-e (OR -Emissivity) 0.85 (OR any value [0, 1]) ]
    + convert .seq to .root 
    + each image is stored in one root file roo/frame_idx.root, idx = [0, ...]
    + an average image is stored in file roo/frame_average.root
    + use optional -e OR -Emissivity to change the emissivity value for the output.

  - ./frameanal.py roo/frame_average.root
    + for help, use: "./frameanal.py -h (OR --help)"
    + the 'config_frame' is automatically generated using configFinder.py
      1. If you do not want to generate a new config file use: ./frameanal.py -mc (OR --manualconfig)
    + 'config_frame' is measured corresponding to your setup. Below are the values configFinder.py finds
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

  - ./defectFinder.py result.root
    + uses cooling pipe result to find flaws on the stave
    + for all fits in the defectFinding algorithm, use ./defectFinder.py -v (OR --verbose) result.root
    + out plots of both pipes with defects found on top in flawplots/
    + out plots fit results (if turned on) for each flaw in flawplots/fits/

* Installation:

  - There are a few softwares needed before running this package. A list below.

  - Python 2.7 and extensions numpy and opencv (needed so that pyROOT will run) 
    + download: python 2.7 can be found many places, it is generally installed
        on mac devices, and can be installed via the 'sudo apt install python'
        command on Ubuntu. numpy and opencv are both extensions and can be installed
        from the same source as python ie 'sudo apt install python-numpy'
    + information: PyRoot is only supported using python 2.7

  - ExifTool (read .seq file produced by Flir IR camera)
    + download: http://www.sno.phy.queensu.ca/~phil/exiftool/
    + information: http://www.sno.phy.queensu.ca/~phil/exiftool/TagNames/FLIR.html (e.g.: .seq to .png)
    + command: "exiftool Test.seq" gives all information about the camera
    + command: "exiftool -Emissivity Test.seq" gives "Emissivity : 0.95"
    + command: "exiftool -FLIR:all Test.seq" gives all properties related to the .seq
    + command: "exiftool -RawThermalImage -b  Test.seq > firstImage.dat" extracting the first image from Test.seq

  - ImageMagick (needed for command 'convert') 
    + download: http://cactuslab.com/imagemagick/

  - ROOT (from CERN to make plots and store data)
    + download: https://root.cern.ch/downloading-root

  - OpenCV (needed for configFinder.py allows for image manipulation)
    + download: https://opencv.org/releases.html or it can also be installed by Homebrew or MacPorts
    + install Homebrew on Mac: /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    + install OpenCV: brew install opencv
      Make sure of doing: mkdir -p /Users/jieyu/Library/Python/2.7/lib/python/site-packages; 
                          echo 'import site; site.addsitedir("/usr/local/lib/python2.7/site-packages")' >> /Users/jieyu/Library/Python/2.7/lib/python/site-packages/homebrew.pth
        after installing opencv


