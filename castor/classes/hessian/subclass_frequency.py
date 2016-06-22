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
    def __init__(self, Thz_in = 0, THz2Pi_in = 0, cm1_in = 0, meV_in = 0, imag_in = False):
        self.THz       = Thz_in
        self.THz2Pi    = THz2Pi_in
        self.cm1       = cm1_in
        self.meV       = meV_in
        self.imaginary = imag_in
        self.atdiff    = dict() # Filled by class_dipol.py during outcar read
        self.intensity = None   # Set by subclass_matrix.py in setIntensities
    def setEigen(self, eigenval_in = 0, atdiff_in = dict(), imag_in = False):
        self.THz       = eigenval_in*15.633304592
        self.THz2Pi    = eigenval_in*98.2269497148
        self.cm1       = eigenval_in*521.47091
        self.meV       = eigenval_in*64.6541499
        self.atdiff    = atdiff_in
        self.imaginary = imag_in
    def getString(self, decimals_in = 6, spaces_in = 3, prefix = None, printmode = "all"):
        self.decimals = decimals_in
        self.spaces   = spaces_in
        string        = ""
        if(prefix!=None):
            string += '{0:>{width}}'.format(prefix, width=self.spaces)
            string += " "
        if(self.imaginary==False):
            string += "f  ="
        elif(self.imaginary==True):
            string += "f/i="
        if(printmode=="all" or printmode=="freq" or self.intensity == None):
            data   = '{0:.{width}f} THz'.format(self.THz, width=self.decimals)
            string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+6)
            data   = '{0:.{width}f} 2PiTHz'.format(self.THz2Pi, width=self.decimals)
            string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+12)
            data   = '{0:.{width}f} cm-1'.format(self.cm1, width=self.decimals)
            string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+11)
            data   = '{0:.{width}f} meV'.format(self.meV, width=self.decimals)
            string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+9)
        if(self.intensity != None and (printmode == "all" or printmode == "intensity")):
            if(printmode == "all"):
                string += "\n"
                data   = '{0:.{width}f} x2'.format(self.intensity[0], width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+13)
                data   = '{0:.{width}f} y2'.format(self.intensity[1], width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+9)
                data   = '{0:.{width}f} z2'.format(self.intensity[2], width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+13)
                data   = '{0:.{width}f} int'.format(self.intensity[3], width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+11)
            else:
                data   = '{0:.{width}f} x2'.format(self.intensity[0], width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+6)
                data   = '{0:.{width}f} y2'.format(self.intensity[1], width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+6)
                data   = '{0:.{width}f} z2'.format(self.intensity[2], width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+6)
                data   = '{0:.{width}f} int'.format(self.intensity[3], width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+7)
        return string
    def write(self, decimals_in = 6, spaces_in = 3, prefix = None, printmode = "all"):
        print self.getString(decimals_in, spaces_in, prefix, printmode)
    def writeDiff(self):
        for diff in self.atdiff:
            print diff