#! /usr/bin/env python

# 
# class_dipol.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2016 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math
import re
import numpy as np

#MY PATHS

#MY CLASSES
from dipols.subclass_matrix import Matrix
from dipols.subclass_frequency import Frequency

class Dipols:
    """Holds the dipols from calculating the Hessian"""
    def __init__(self):
        self.decimals  = 6
        self.spaces    = 3
        self.matrix    = None
        self.finitdiff = None
        self.elnr      = np.dtype(int)    # Vector of element numbers.
        self.elements  = []               # List of element names.
        self.masses    = []               # List of element masses.
    def read(self, outcar_in):
        readmode = 0
        with open(outcar_in, 'r') as inputfile:
            self.nions           = None
            self.dipol_diff      = []
            self.frequencies     = []
            self.degrees_freedom = []
            vaspversion          = None
            nwrite               = None   # Which nwrite has been used.
            dipol_correction     = None
            self.dipol_store     = None
            self.dipol_ref       = None
            scale_diff           = False  # Whether displacements should be corrected for mass.
            for line in inputfile:
                if(readmode == 0):
                    #which version of vasp has been used?
                    match = re.search('^[ \t]*vasp[.]([0-9])[.].*$',line)
                    if(match): # vasp version found
                        vaspversion = int(match.group(1))
                        if(vaspversion != 4 and vaspversion != 5):
                            print "Error: Unknown Vasp version: " + str(vaspversion)
                            sys.exit()
                        readmode += 1
                        continue
                elif(readmode == 1): # Search for elemental properties
                    # Which element name?
                    match = re.search('^[ \t]*TITEL[ \t]+=[ \t]+[A-Z_]+ ([A-Za-z]+) [0-9A-Za-z]+[ \t]*$',line)
                    if(match):
                        self.elements.append(match.group(1))
                        continue
                    # What element mass?
                    match = re.search('^[ \t]*POMASS[ \t]+=[ \t]+([0-9.]+);.*$',line)
                    if(match):
                        self.masses.append(float(match.group(1)))
                        continue
                    # How many atoms per type of element?
                    match = re.search('^[ \t]*ions per type =([0-9 \t]+)[ \t]*$',line)
                    if(match):
                        stringlist_values = match.group(1).split()
                        intlist_values = []
                        for i in range(0,len(stringlist_values)):
                            intlist_values.append(int(stringlist_values[i]))
                        self.elnr = np.array(intlist_values)
                        continue
                    # Which NWRITE is used?
                    match = re.search('^[ \t]*NWRITE[ \t]*=[ \t]*([0-4])[ \t]+.*$',line)
                    if(match):
                        nwrite   = int(match.group(1))
                        readmode += 1
                        continue
                elif(readmode == 2):
                    #which LDIPOL tag has been set
                    match = re.search('^[ \t]*LDIPOL[ \t]*=[ \t]*([FT])[ \t]*.*$',line)
                    if(match): # LDIPOL tag found
                        if(match.group(1)=="T"):
                            dipol_correction = True
                        elif(match.group(1)=="F"):
                            dipol_correction = False
                        else:
                            print "Error: Unknown LDIPOL tag: " + line
                            sys.exit()
                        readmode += 1
                        continue
                elif(readmode == 3):
                    #read dipole moments
                    match = re.search('^[ \t]*dipolmoment[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*electrons.*$',line)
                    # " dipolmoment          -0.000002     -0.051875     -0.352631 electrons x Angstroem"
                    if(match): # dipol moment found
                        self.dipol_store = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                        if(dipol_correction == False):
                            if(self.dipol_ref==None):
                                self.dipol_ref = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                            else:
                                dipol = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                                self.dipol_diff.append(dipol)
                        continue

                    #for LDIPOL = T, dipols are found for every electronic step, so only use the converged step
                    if(dipol_correction == True):
                        match = re.search('^.*aborting loop because EDIFF is reached.*$',line)
                        # "------------------------ aborting loop because EDIFF is reached ----------------------------------------"
                        if(match): # electronic convergence
                            if(self.dipol_ref==None):
                                self.dipol_ref   = self.dipol_store
                                self.dipol_store = None
                            else:
                                self.dipol_diff.append(self.dipol_store)
                                self.dipol_store = None
                            continue
                
                    #read the degrees of freedom
                    # "              11X         11Y         11Z         12X         12Y         12Z"
                    match = re.search('^([ \t]+[0-9]+[XYZ])+[ \t]*$',line)
                    if(match):
                        self.degrees_freedom = match.group(0).split()
                        readmode += 1
                        continue
                elif(readmode == 4 or readmode == 5):
                    #make sure we are reading the correct dynamical matrix
                    if(re.search('^[ \t]*Eigenvectors and eigenvalues of the dynamical matrix[ \t]*$',line)):
                        readmode += 1
                        if(vaspversion == 5 and nwrite < 3):
                            scale_diff = True
                            readmode   += 1
                        continue
                elif(readmode == 6):
                    match = re.search('^[ \t]*[0-9]+ f  =[ \t]*([0-9.]+) THz[ \t]*([0-9.]+) 2PiTHz[ \t]*([0-9.]+) cm-1[ \t]*([0-9.]+) meV[ \t]*$',line)
                    # "6 f  =    2.865417 THz    18.003946 2PiTHz   95.580018 cm-1    11.850416 meV"
                    if(match): # real freq found
                        atnum = 0
                        self.frequencies.append(Frequency(float(match.group(1)), float(match.group(2)), float(match.group(3)), float(match.group(4)),False))
                        continue
                    match = re.search('^[ \t]*[0-9]+ f/i=[ \t]*([0-9.]+) THz[ \t]*([0-9.]+) 2PiTHz[ \t]*([0-9.]+) cm-1[ \t]*([0-9.]+) meV[ \t]*$',line)
                    if(match): # imag freq found
                        atnum = 0
                        self.frequencies.append(Frequency(float(match.group(1)), float(match.group(2)), float(match.group(3)), float(match.group(4)),True))
                        continue
                    if(len(self.frequencies)>0):
                        match = re.search('^[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*$',line)
                        if(match): # reading atomic displacement line
                            atnum += 1
                            diff = [float(match.group(4)), float(match.group(5)), float(match.group(6))]
                            if(scale_diff == True):
                                mass    = self.getMass(atnum)
                                diff[0] /= math.sqrt(mass)
                                diff[1] /= math.sqrt(mass)
                                diff[2] /= math.sqrt(mass)
                            self.frequencies[-1].atdiff.append(diff)
                            continue
                        match = re.search('^[ \t]*Finite differences POTIM=[ \t]*([0-9.E-]+)[ \t]*$',line)
                        # "Finite differences POTIM=  2.00000000000000004E-002"
                        if(match): # reading POTIM and ending frequency reading
                            self.finitdiff = float(match.group(1))
                            readmode += 1 # making sure no more freqs are read
                            continue
                elif(readmode == 7):
                    break

        if(readmode == 0):
            print "Error: Could not find Vasp version"
        elif(readmode ==1):
            print "Error: Could not find elemental properties"
        elif(readmode ==2):
            print "Error: Could not find LDIPOL tag"
        elif(readmode ==3):
            print "Error: Could not find dipol moments"
        elif(readmode ==4 or readmode == 5):
            print "Error: Could not find dynamical matrix"
        elif(readmode ==6):
            print "Error: Could not find all frequencies"
    def getMass(self, number_in):
        atomsum = 0
        for number, mass in zip(self.elnr, self.masses):
            atomsum += number
            if(number_in <= atomsum):
                return mass
        print "Could not get element mass"
        sys.exit()
    def getMatrix(self):
        self.matrix = Matrix()
        self.matrix.setup(self.dipol_diff,self.degrees_freedom,self.finitdiff)
        self.matrix.setIntensities(self.frequencies)
    def writeList3(self, list3, decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces   = spaces_in
        string        = ""
        xcor          = '{0:.{width}f}'.format(list3[0], width=self.decimals)
        ycor          = '{0:.{width}f}'.format(list3[1], width=self.decimals)
        zcor          = '{0:.{width}f}'.format(list3[2], width=self.decimals)
        string        += '{0:>{width}}'.format(xcor, width=self.decimals+self.spaces+4)
        string        += '{0:>{width}}'.format(ycor, width=self.decimals+self.spaces+4)
        string        += '{0:>{width}}'.format(zcor, width=self.decimals+self.spaces+4)
        return string
    def write(self,printmode = "all"):
        if(printmode=="all"):
            print "Initial dipole moment: " + self.writeList3(self.dipol_ref)
            self.matrix.printMatrix()
        count = 0
        for freq in self.frequencies:
            count += 1
            freq.write(prefix=count, printmode=printmode)