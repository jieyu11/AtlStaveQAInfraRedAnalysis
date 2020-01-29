#!/usr/bin/env python
import cv2
import numpy as np
from matplotlib import pyplot as plt



def filter(x,cut):
    if x < cut:
        return 0
    else:
        return x



img = 10*cv2.imread('image.tiff',0)

np.apply_over_axis(lambda x : x filter(x,58), 0, img)
np.app


'''
kernel = np.ones((5,5),np.float32)/25
img = cv2.filter2D(img,-1,kernel)
'''

edges = cv2.Canny(img,0,150)
#edges = cv2.Canny(img_smooth,50,150)


minLineLength = 100
maxLineGap = 20
#lines = cv2.HoughLinesP(edges,5,np.pi/180,100,minLineLength,maxLineGap)
# vertical lines:
lines = cv2.HoughLinesP(edges,rho = 1,theta = 1*np.pi/10000,threshold = 20,minLineLength = 10, maxLineGap = 5)
#horizontal lines:
#lines = cv2.HoughLinesP(edges,rho = 1,theta = 1*np.pi/1000,threshold = 100,minLineLength = 200,maxLineGap = 175)

for line in lines:
    for x1,y1,x2,y2 in line:
        cv2.line(img,(x1,y1),(x2,y2),(255,0,0),2)


plt.subplot(121),plt.imshow(img)
plt.subplot(122),plt.imshow(edges,cmap = 'gray')
plt.title('Edge Image'), plt.xticks([]), plt.yticks([])

plt.show()
