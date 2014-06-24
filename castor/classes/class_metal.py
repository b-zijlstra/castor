#! /usr/bin/env python

# 
# class_metal.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys

#MY PATHS

#MY CLASSES

class Metal:
    """Defines a metal packing"""
    def __init__(self, metal_in, lc_in, hc_in, packing_in):
        self.metal = metal_in
        self.lc = float(lc_in)
        self.hc = float(hc_in)
        self.packing = packing_in