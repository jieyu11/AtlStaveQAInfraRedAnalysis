#!/usr/bin/env python
import configparser

config = configparser.ConfigParser()

config.read("parameters.cfg")

Tin = config["Default"]["Tin"]
Tout = config["Default"]["Tout"]
Cliq = config["Default"]["Cliq"]
Zcut = config["Default"]["Zcut"]
FR = config["Default"]["FR"]

#print(config.items("Default"))

for variable in config.items("Default"):
  print(variable[0] + " = " + variable[1]) 

