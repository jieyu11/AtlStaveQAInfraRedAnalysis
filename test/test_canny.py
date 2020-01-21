#!/usr/bin/env python
import cv2
import numpy as np
from matplotlib import pyplot as plt

img = 10*cv2.imread('image.tiff',0)


edges = cv2.Canny(img,50,150)

minLineLength = 100
maxLineGap = 100
lines = cv2.HoughLinesP(edges,1,np.pi/180,100,minLineLength,maxLineGap)

for x1,y1,x2,y2 in lines[0]:
    cv2.line(img,(x1,y1),(x2,y2),(255,0,0),2)


plt.subplot(121),plt.imshow(img)
plt.title('Original Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(edges,cmap = 'gray')
plt.title('Edge Image'), plt.xticks([]), plt.yticks([])

plt.show()
