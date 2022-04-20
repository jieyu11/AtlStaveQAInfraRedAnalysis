# ITk Stave Impedance Quality Control

This code takes an IR image (in the CSV format) of a stave and outputs the thermal impedances of the thermal heat paths between the cooling pipe and the surface in selected regions.


## Dependencies
* Python >= 3.8

## Installation

First initialize the virtual environment

```bash
python -m venv ./env
```
and activate it
```bash
source env/bin/activate
```
Install the required modules
```bash
pip install -r requirements.txt
```
The virtual environment can be deactivated via
```bash
deactivate
```

## Usage
Don't forget to activate the virtual environment before running the script
```bash
source env/bin/activate
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
