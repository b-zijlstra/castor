#! /usr/bin/env python

# 
# subclass_frequency.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math

#MY PATHS

#MY CLASSES

class Frequency:
    """Defines a frequency"""
    def __init__(self, eigenval_in = 0, imag_in = False):
        self.THz = eigenval_in*15.633304592
        self.THz2Pi = eigenval_in*98.2269497148
        self.cm1 = eigenval_in*521.47091
        self.meV = eigenval_in*64.6541499
        self.imaginary = imag_in
    def getString(self, decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces = spaces_in
        string = ""
        if(self.imaginary==False):
            string += "f  ="
        elif(self.imaginary==True):
            string += "f/i="
        data = '{:.{width}f} THz'.format(self.THz, width=self.decimals)
        string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+6)
        data = '{:.{width}f} 2PiTHz'.format(self.THz2Pi, width=self.decimals)
        string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+12)
        data = '{:.{width}f} cm-1'.format(self.cm1, width=self.decimals)
        string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+11)
        data = '{:.{width}f} meV'.format(self.meV, width=self.decimals)
        string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+9)
        return string
    def write(self, decimals_in = 6, spaces_in = 3):
        print self.getString(decimals_in, spaces_in)