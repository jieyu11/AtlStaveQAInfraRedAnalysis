#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 15:18:14 2020

@author: vozdecky
"""

import cv2
from matplotlib import pyplot as plt
import numpy as np
import csv

image = []
with open('CSV_with_0.92_emmy.csv') as csvfile:
  reader = csv.reader(csvfile, quoting=csv.QUOTE_NONNUMERIC)
  for row in reader:
    image.append(row)

img = np.array(image)

height = img.shape[0]
width = img.shape[1]

upper_img = img[0:int(height/2), 0:width].transpose()
lower_img = img[int(height/2):height, 0:width].transpose()



test = np.ones((600,400))*255

img2 = (img - np.min(img))/(np.max(img)-np.min(img))*200 + 50 
cv2.rectangle(img2,(100,100),(200,200),(0,0,0),thickness=3)

plt.imshow(img2)
plt.show()
