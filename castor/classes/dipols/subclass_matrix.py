#! /usr/bin/env python

# 
# subclass_matrix.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2016 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math
import numpy as np

#MY PATHS

#MY CLASSES
from subclass_frequency import Frequency

class Matrix:
    """Defines the dipol matrix from a Hessian"""
    def __init__(self):
        self.decimals = 6
        self.spaces   = 3
        self.dipols   = []
    def setup(self, dipol_diff, degrees_freedom, diff):
        #check if the amount of dipols match the degrees of freedom
        if(len(dipol_diff)!= 2*len(degrees_freedom)):
            print "Error, amounts of dipols and degrees of freedom do not match!"
            sys.exit()

        atom      = 0
        direction = 0
        dipcount  = 0
        for degree in degrees_freedom:
            atom      = int(degree[:-1])
            direction = degree[-1]
            #first dipol_diff is +diff, second dipol_diff is -diff
            mux = (dipol_diff[dipcount][0]-dipol_diff[dipcount+1][0])/(2*diff)
            muy = (dipol_diff[dipcount][1]-dipol_diff[dipcount+1][1])/(2*diff)
            muz = (dipol_diff[dipcount][2]-dipol_diff[dipcount+1][2])/(2*diff)
            self.dipols.append([atom,direction,mux, muy, muz])
            dipcount += 2
    def setIntensities(self,frequencies):
        for freq in frequencies:
            sumx    = 0
            sumy    = 0
            sumz    = 0
            atcount = 0
            for atom in freq.atdiff:
                atcount += 1
                for dip in self.dipols:
                    if(dip[0] == atcount):
                        if(dip[1] == "X"):
                            sumx+=dip[2]* atom[0]
                            sumy+=dip[3]* atom[0]
                            sumz+=dip[4]* atom[0]
                        elif(dip[1] == "Y"):
                            sumx+=dip[2]* atom[1]
                            sumy+=dip[3]* atom[1]
                            sumz+=dip[4]* atom[1]
                        elif(dip[1] == "Z"):
                            sumx+=dip[2]* atom[2]
                            sumy+=dip[3]* atom[2]
                            sumz+=dip[4]* atom[2]
            totaldipol = sumx**2 + sumy**2 + sumz**2
            freq.intensity = [sumx**2,sumy**2,sumz**2,totaldipol]

    def printMatrix(self):
        for dipol in self.dipols:
            print dipol

    def printMatrix(self, decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces = spaces_in
        for dipol in self.dipols:
            string = ""
            string += '{0:>{width}}'.format(dipol[0], width=self.spaces)
            string += '{0:>{width}}'.format(dipol[1], width=self.spaces)
            xcor   = '{0:.{width}f}'.format(dipol[2], width=self.decimals)
            ycor   = '{0:.{width}f}'.format(dipol[3], width=self.decimals)
            zcor   = '{0:.{width}f}'.format(dipol[4], width=self.decimals)
            string += '{0:>{width}}'.format(xcor, width=self.decimals+self.spaces+4)
            string += '{0:>{width}}'.format(ycor, width=self.decimals+self.spaces+4)
            string += '{0:>{width}}'.format(zcor, width=self.decimals+self.spaces+4)
            print string