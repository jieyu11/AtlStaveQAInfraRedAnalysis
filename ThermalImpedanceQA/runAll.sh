#!/bin/bash

# This script runs the python script in impedancesFromCSV.py on every .csv file in the /data folder

if [ ! -d ./output ]
    then
    echo "The output folder doesn't exist. Creating /output."
    mkdir output
fi

for file in `ls data/*.csv`; do
    ./impedanceFromCSV.py $file
done

