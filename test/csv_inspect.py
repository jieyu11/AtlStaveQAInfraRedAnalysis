#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:40:33 2020

@author: vozdecky
"""

import csv
import numpy as np
from matplotlib import pyplot as plt

results = []
with open("CSV_with_0.92_emmy.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        results.append(row)

results2 = []
with open("image_from_raw.csv") as csvfile:
    reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC) # change contents to floats
    for row in reader: # each row is a list
        results2.append(row)

diff = abs(np.array(results) - np.array(results2).transpose())

plt.imshow(diff, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()
