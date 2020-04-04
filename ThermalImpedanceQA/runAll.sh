#!/bin/bash

# This script runs the python script in impedancesFromCSV.py on every .csv file in the /data folder

for file in `ls data/*.csv`; do
    ./impedanceFromCSV.py $file
done

