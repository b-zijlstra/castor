#! /usr/bin/env python

# 
# class_hessian.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math
import re
import numpy as np

#MY PATHS

#MY CLASSES

class Chargemol:
    """Defines chargemol output data"""
    def __init__(self):
        # Storage
        self.BOsums       = []   # sums of bond orders
        self.pairmatrix   = []   # matrix of bond pairs and related data

    def read(self, input_in):
        with open(input_in, 'r') as inputfile:
            readmode = 0
            for line in inputfile: # Read the chargemol output line for line
                if(readmode == 0): # Check for sum of BOs
                    match = re.search('^[ \t]* The sum of BOs for each atom are:[ \t]*$',line)
                    if(match): # sums found
                        print "FOUND BO SUMS"
                        readmode += 1
                        continue
                elif(readmode == 1): # placeholder
                    # placeholder
                    match = re.search('^[ \t]*([0-9.E-]+)[ \t]*$',line)
                    if(match):
                        self.BOsums.append(float(match.group(1)))
                        continue
                    match = re.search('^[ \t]*$',line)
                    if(match): # end of sums
                        placeholder     = None # placeholder
                        readmode += 1
                        continue
                elif(readmode == 2): # placeholder
                    # placeholder
                    match = re.search('^[ \t]*TITEL[ \t]+=[ \t]+[A-Z_]+ ([A-Za-z]+) [0-9A-Za-z]+[ \t]*$',line)
                    if(match):
                        placeholder     = None # placeholder
                        continue
                    match = re.search('^[ \t]*POTIM[ \t]*=[ \t]*([0-9.E-]+)[ \t]+.*$',line)
                    if(match):
                        placeholder     = None # placeholder
                        readmode += 1
                        continue
                elif(readmode == 3): # placeholder
                    # placeholder
                    match = re.search('^[ \t]*TITEL[ \t]+=[ \t]+[A-Z_]+ ([A-Za-z]+) [0-9A-Za-z]+[ \t]*$',line)
                    if(match):
                        placeholder     = None # placeholder
                        continue
                    match = re.search('^[ \t]*POTIM[ \t]*=[ \t]*([0-9.E-]+)[ \t]+.*$',line)
                    if(match):
                        placeholder     = None # placeholder
                        readmode += 1
                        continue
                elif(readmode == 4):
                    break
        if(readmode == 0):
            print "Error: 0"
        elif(readmode ==1):
            print "Error: 1"
        elif(readmode ==2):
            print "Error: 2"
        elif(readmode ==3):
            print "Error: 3"
    def getbond(self, atoms1_in, atoms2_in, printmode):

        if(printmode == "normal" or printmode == "all"):
            for atom in atoms1_in:
                print "BOsum of " + str(atom) + " = " + str(self.BOsums[atom]) 
            if(printmode == "all"):
                print "BO sums:"
                for BOsum in self.BOsums:
                    print BOsum
        else:
            placeholder     = None # placeholder
            if(printmode == "least"):
                placeholder     = None # placeholder
            else:
                placeholder     = None # placeholder
        print "Hello World!"