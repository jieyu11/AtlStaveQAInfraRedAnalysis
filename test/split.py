#!/usr/bin/env python
import sys
import cv2

#crops an image and leaves the upper half

path = sys.argv[1]

img = cv2.imread(path)

height = img.shape[0]
width = img.shape[1]

lower_img = img[0:int(height/2), 0:width]
upper_img = img[int(height/2):height, 0:width]

cv2.imwrite(path, upper_img)

print("SPLITING!!!!!")


'''
height = img.shape[0]
width = img.shape[1]

lower_img = img[0:int(height/2), 0:width]
upper_img = img[int(height/2):height, 0:width]

cv2.imwrite('pic_lower.jpg', upper_img)
cv2.imwrite('pic_upper.jpg', lower_img)
'''
