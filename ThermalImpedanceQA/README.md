# ITk Stave Impedance Quality Control

This code takes an IR image (in the CSV format) of a stave and outputs the thermal impedances of the thermal heat paths between the cooling pipe and the surface in selected regions.


## Dependencies
* Anaconda 3 (for pytest, numpy, configparsers etc.)
* OpenCV

## Installation

Download and install [Anaconda](https://www.anaconda.com/products/individual) to manage the dependencies. Make sure Anaconda is added to the system PATH variable.

First initialize Anaconda:

```bash
conda init bash
```
And set the auto-activation to false to prevent Anaconda initializing automatically whenever a new shell window is opened:
```bash
conda config --set auto_activate_base false
```
Restart the shell.
The evironment can be now activated through:
```bash
conda actiavte
```

You can check the list of all modules installed
```bash
conda list
```
and make sure they are all up-to-date with
```bash
conda update --all
```
Some required packages won't be in the default Anaconda envirnoment and they need to be installed:
```bash
conda install opencv
```
Activate the executables:
```bash
chmod +x *.py
```
Finally make sure your python version is 3.8.5 or newer
```bash
python --version
```
and all dependencies are correclty installed:
```bash
./test_dependencies.py
```
If you get no errors you are good to go.

At the end you may want to deactivate Anaconda:
```bash
conda deactivate
```

## Usage
After opening a shell make sure you activate the envirnoment
```bash
conda actiavte
```

When you are done you can deactivate the environment
```bash
conda deactivate
```

```
usage: ./impedanceFromCSV.py [-h] [-g] [-1f] [-d]
                          path_to_image parameters_file
```

Required arguments:

| __Argument__  | __Description__                                                      |
| :---          | :---                                                                 |
| `path_to_image` | path to the csv file containing the thermal image of the stave |
| `parameters_file` | path to the .cfg file with the configuration parameters |

Optional arguments:

| __Argument__         |  __Description__                                                                                                                                   |   
| :---                 | :---                                                                                                                                               |   
| `-h, --help`         | Show the help message and exit.                                                                                                                    |   
| `-d, --debug`        | Run the script in the debug mode. It will output a debug log into the debug_output subdirectory.                                                   |   
| `-g --graphs`        | Produce graphs.                                                                                                                                    |
| `-1f, --one_face`    | Use this flag when there is only one face in the IR image. By default the code will look for two faces produced by the Queen Mary V-mirror set-up. |


Examples:
```bash
./impedanceFromCSV.py data/CSV_with_0.92_emmy.csv parameters/water-default.cfg -g
```
```bash
./impedanceFromCSV.py data/FEA_with_honeycomb.csv parameters/water-default.cfg -1f -g 
```

## Testing
Before committing any changes make sure the modified version passes all the tests by runnning pytest in this folder:
```bash
pytest
```

## Contact
If you have any questions, please do not hesitate to contact me via email at l.vozdecky@qmul.ac.uk.
