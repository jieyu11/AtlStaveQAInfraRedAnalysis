# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 09:30:28 2020

Testing script for the Stave class.

@author: LuBoss
"""
import pytest
import numpy as np
from stave import Region, GeneralRegion

def test_region():
  global_image = 18*np.ones([100,200])
  region_image = 30*np.ones([50,50])
  global_image[25:75,100:150] = region_image
  
  region_30 = Region(global_image,110,140,30,70)
  region_18 = Region(global_image,10,60,20,75)
  
  assert region_30.getAverageTemperature() == 30
  assert region_18.getAverageTemperature() == 18

  assert region_30.getPosition() == [110,140,30,70]
  assert region_18.getPosition() == [10,60,20,75]
  

def test_generalRegion():
  global_image = 18*np.ones([100,200])
  region_image = 30*np.ones([50,50])
  global_image[25:75,100:150] = region_image
  
  region_definition_30 = np.zeros([100,200])
  region_definition_30[30:70,110:140] = np.ones([40,30])
  region_definition_18 = np.zeros([100,200])
  region_definition_18[10:60,20:75] = np.ones([50,55])
  
  
  region_30 = GeneralRegion(global_image,region_definition_30)
  region_18 = GeneralRegion(global_image,region_definition_18)
  
  assert region_30.getAverageTemperature() == 30
  assert region_18.getAverageTemperature() == 18

"""
@pytest.fixture
def wrongRatioStave():
  global_image = 18*np.ones([1240,3160])
  region_image = 30*np.ones([270,3000])
  global_image[400:670,80:3080] = region_image
  return global_image

@pytest.fixture
def correctRatioStave():
  global_image = 18*np.ones([1240,3160])
  region_image = 30*np.ones([251,3000])
  global_image[400:651,80:3080] = region_image
  mystave = Stave(global_image, "parameters.cfg")
  return global_image

def test_ratio(wrongRatioStave,correctRatioStave):
"""