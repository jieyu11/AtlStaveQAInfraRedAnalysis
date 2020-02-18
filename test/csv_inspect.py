#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:40:33 2020

@author: vozdecky
"""

import csv
import numpy as np
import cv2
from matplotlib import pyplot as plt


results = []
with open("A655-000277-312_16_15_05_892.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        results.append(row)

results2 = []
with open("A655-000277-312_16_15_05_892_RAW.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        results2.append(row)


#images must be shifted as they are offset
diff = abs(np.array(results) - np.rot90(np.array(results2)))
ratio = np.array(results) / np.rot90(np.array(results2))
plt.imshow(diff, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()
