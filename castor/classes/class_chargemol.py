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
                        readmode += 1
                        continue
                elif(readmode == 1): # Read the BO sums
                    match = re.search('^[ \t]*([0-9.E-]+)[ \t]*$',line)
                    if(match):
                        self.BOsums.append(float(match.group(1)))
                        continue
                    match = re.search('^[ \t]*$',line)
                    if(match): # end of sums
                        readmode += 1
                        continue
                elif(readmode == 2): # Check for start of pair matrix
                    match = re.search('^[ \t]* The final bond pair matrix is*$',line)
                    if(match):
                        readmode += 1
                        continue
                elif(readmode == 3): # Read the pair matrix
                    match = re.search('^[ \t]* The legend for the bond pair matrix follows:.*$',line)
                    if(match):
                        readmode += 1
                        continue
                    else:
                        numberstring = line
                        numberstring = numberstring.split()
                        row = []
                        try:
                            counter = 0
                            for number in numberstring:
                                counter += 1
                                if counter <= 2:
                                    row.append(int(number))
                                else:
                                    row.append(float(number))
                        except ValueError:
                            print "ValueError - invalid pair matrix line:" + line
                            self.help()
                            sys.exit()
                        self.pairmatrix.append(row)
                        continue
                elif(readmode == 4):
                    break
        if(readmode == 0):
            print "Error: Could not find bond order sums"
        elif(readmode ==1):
            print "Error: Could not read bond order sums"
        elif(readmode ==2):
            print "Error: Could not find pair matrix"
        elif(readmode ==3):
            print "Error: Could not read pair matrix"
    def getbond(self, atoms1_in, atoms2_in, printmode):
        if(printmode == "normal" or printmode == "all"):
            for atom in atoms1_in:
                BOsum = self.BOsums[atom]
                BOsum_reduced = BOsum
                print "BOsum of " + str(atom) + " = " + str(BOsum) 
                for otheratom in atoms2_in:
                    BO = 0.0
                    for row in self.pairmatrix:
                        if (row[0] == atom and row[1] == otheratom) or (row[0] == otheratom and row[1] == atom):
                            if row[19] > BO:
                                BO = row[19]
                    BOsum_reduced -= BO
                    print "BO of " + str(atom) + " with " + str(otheratom) + " = " + str(BO)
                print "BOsum - BO others = " + str(BOsum_reduced) 
            if(printmode == "all"):
                print "BO sums:"
                for BOsum in self.BOsums:
                    print BOsum
                print "pair matrix:"
                for row in self.pairmatrix:
                    rowstring = ""
                    rowstring += '{0:>{width}}'.format(row[0], width=5)
                    rowstring += '{0:>{width}}'.format(row[1], width=5)
                    rowstring += '{0:>{width}}'.format(row[2], width=7)
                    rowstring += '{0:>{width}}'.format(row[3], width=7)
                    rowstring += '{0:>{width}}'.format(row[4], width=7)
                    rowstring += '{0:>{width}}'.format(row[19], width=12)
                    print rowstring
        else:
            placeholder     = None # placeholder
            if(printmode == "least"):
                placeholder     = None # placeholder
            else:
                placeholder     = None # placeholder