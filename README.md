* Purpose:
  
  - Convert raw ADC counts recorded by the Flir InfraRed Camera into temperature (C).
  - Analyze the stave temperature profile and extract the temperature profile of the cooling pipe in the stave.

* Run:
  - Check installation below first!

  - Once you believe that the system is set appropriately try running on the example file

    ./seqToProfile.py FlirIR_rec_example.seq

    This will run the read_sequence.sh and frameanal.py scripts for the FlirIR_rec_example.seq file.
    This should generate a folder called: plot-FlirIR_rec_example/   containing 
    two more folders called: hist/ and fit/  along with plots of the stave core.

    Next test the defect finder.

    ./defectFinder.py plot-FlirIR_rec_example/result.root

    This python script will plot all attached result.root files and can be used
    to make all sorts of plots.

  - ./seqToProfile.py [sequence files]
    + performs the full conversion on each applied sequence file. If it fails to
      find a configuration it may not save the results and you may have to individually
      do the steps below.

      1. performs read_sequence with preset emissivity 0.92
      2. performs frameanal.py
      3. reorganizes all of the appropriate needed files into a single folder
         based off the name of each individual sequence file.

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
    (Most(I think all) of these packages can be easily downloaded in Ubuntu using apt-get install [package])

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

+ Other Packages

  - SliceFinder.py
    + This package takes a stave core's average_frame.root file and then converts
      it into average strips across each of the modules. 
    + It is currently written for a 13 module stave core.
    + To Run with example files:
      ./extras/SliceFinder.py extras/SliceFinder_Example_Files/Example_Stave7L_m55/frame_average.root
    + This will use an example .csv file from graham with 14 modules(it ignores 
      the ones closest to the end of stave card (I Think).
    + This also requires scipy to work. It is another python package that you may not be standard.

+ Questions?
  - Send an email to wheidorn@iastate.edu. I may be able to help.

