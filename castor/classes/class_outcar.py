#! /usr/bin/env python3

# 
# class_outcar.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2017 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math
import re
import numpy as np

#MY PATHS

#MY CLASSES
#from class_vector import Vector
from .class_unitcell import Unitcell

class Outcar:
    """Defines an OUTCAR"""
    def __init__(self):
        # General storage
        self.ibrion      = None             # IBRION tag.
        self.decimals    = 6                # Number of decimals to output.
        self.spaces      = 3                # Number of spaces to output.
        self.elnr        = np.dtype(int)    # Vector of element numbers.
        self.elements    = []               # List of element names.
        self.finitdiff   = None             # Finite difference used.

        # Storage for Dipol matrix
        self.volume      = None             # Volume of the cell.
        self.idipol      = None             # Which IDIPOL has been used.
        self.ldipol      = None             # Whether dipol corrections are included in the potential.
        self.epsilon     = None             # Bulk dielectric constant.
                                            # 
        # Storage for the frames
        self.numframes   = 0                # Number of frames.
        self.positions   = []               # List of positions for each frame.
        self.forces      = []               # List of forces for each frame.
        self.dipols      = []               # List of dipols for each frame.
        self.cell        = []               # List of unitcell vectors for each frame.

    def read(self, outcar_in):
        readmode = 0
        sublevel = None
        with open(outcar_in, 'r') as inputfile:
            # OUTCAR info
            vaspversion = None   # Which Vasp version has been used.
            nwrite      = None   # Which nwrite has been used.
            self.idipol = None   # Which IDIPOL has been used.
            self.ldipol = None   # Whether dipol corrections are included in the potential.

            # Temporary storage for dipols
            dipol_store     = None # Storage for the last found dipol moment.

            for line in inputfile: # Read the OUTCAR line for line
                if(readmode == 0): # Get OUTCAR info
                    # Which version of vasp has been used?
                    match = re.search('^[ \t]*vasp[.]([0-9])[.].*$',line)
                    if(match): # vasp version found
                        vaspversion = int(match.group(1))
                        if(vaspversion != 4 and vaspversion != 5):
                            print("Error: Unknown Vasp version: " + str(vaspversion))
                            sys.exit()
                        readmode += 1
                        continue
                elif(readmode == 1): # Search for elemental properties
                    # Which element name?
                    match = re.search('^[ \t]*TITEL[ \t]+=[ \t]+[A-Z_]+ ([A-Za-z]+) [0-9A-Za-z]+[ \t]*$',line)
                    if(match):
                        self.elements.append(match.group(1))
                        continue
                    # How many atoms per type of element?
                    match = re.search('^[ \t]*ions per type =([0-9 \t]+)[ \t]*$',line)
                    if(match):
                        stringlist_values = match.group(1).split()
                        intlist_values    = []
                        for i in range(0,len(stringlist_values)):
                            intlist_values.append(int(stringlist_values[i]))
                        self.elnr = np.array(intlist_values)
                        continue
                    # Which NWRITE is used?
                    match = re.search('^[ \t]*NWRITE[ \t]*=[ \t]*([0-4])[ \t]+.*$',line)
                    if(match):
                        nwrite   = int(match.group(1))
                        continue
                    # What IBRION was used?
                    match = re.search('^[ \t]*IBRION[ \t]*=[ \t]*([0-9])[ \t]+.*$',line)
                    if(match):
                        self.ibrion = int(match.group(1))
                        continue
                    # What POTIM was used?
                    match = re.search('^[ \t]*POTIM[ \t]*=[ \t]*([0-9.E-]+)[ \t]+.*$',line)
                    if(match):
                        self.finitdiff = float(match.group(1))
                        readmode += 1
                        continue
                elif(readmode == 2): # Read the DIPOL tags.
                    # Which IDIPOL tag has been set
                    if(self.idipol == None):
                        match = re.search('^[ \t]*IDIPOL[ \t]*=[ \t]*([0-4])[ \t]+.*$',line)
                        if(match): # IDIPOL tag found
                            if(int(match.group(1)) >= 0 and int(match.group(1)) <= 4):
                                self.idipol = int(match.group(1))
                            else:
                                print("Error: Unknown IDIPOL tag: " + line)
                                sys.exit()
                            continue
                    # Which LDIPOL tag has been set
                    if(self.ldipol == None):
                        match = re.search('^[ \t]*LDIPOL[ \t]*=[ \t]*([FT])[ \t]+.*$',line)
                        if(match): # LDIPOL tag found
                            if(match.group(1)=="T"):
                                self.ldipol = True
                            elif(match.group(1)=="F"):
                                self.ldipol = False
                            else:
                                print("Error: Unknown LDIPOL tag: " + line)
                                sys.exit()
                            continue
                    # What EPSILON has been set
                    if(self.epsilon == None):
                        match = re.search('^[ \t]*EPSILON[ \t]*=[ \t]*([0-9.E]+)[ \t]+.*$',line)
                        if(match): # EPSIOLON tag found
                            self.epsilon = float(match.group(1))
                            continue
                    # What is the volume of the cell
                    if(self.volume == None):
                        match = re.search('^[ \t]*volume of cell[ \t]*:[ \t]*([0-9.]+)[ \t]*$',line)
                        if(match): # volume found
                            self.volume = float(match.group(1))
                            if(self.epsilon == None):
                                self.epsilon = 1.0
                            readmode += 1
                            continue
                elif(readmode == 3): # Read the positions and forces and read the dipolmoments if IDIPOL > 0
                    if(sublevel == "cell"):
                        match = re.search('^.*length of vectors.*$',line)
                        # "  length of vectors"
                        if(match): # all positions and forces read
                            sublevel = None
                            self.cell.append(cell_frame)
                            continue
                        match = re.search('^[ \t]*([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]*$',line)
                        # "     7.687311615 -3.937386022  0.000000000     0.130084489  0.000000000  0.000000000"
                        if(match): # position and force found
                            cell_frame.append(np.array((float(match.group(1)), float(match.group(2)), float(match.group(3)))))
                    elif(sublevel == "forces"):
                        match = re.search('^.*total drift:.*$',line)
                        # "    total drift:                               -0.002172     -0.000568      0.002532"
                        if(match): # all positions and forces read
                            sublevel = None
                            self.numframes += 1
                            self.positions.append(positions_frame)
                            self.forces.append(forces_frame)
                            continue
                        match = re.search('^[ \t]*([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]*$',line)
                        # "      0.72358      1.25327      9.49588         0.001132      0.001590      0.000119"
                        if(match): # position and force found
                            positions_frame.append(np.array((float(match.group(1)), float(match.group(2)), float(match.group(3)))))
                            forces_frame.append(np.array((float(match.group(4)), float(match.group(5)), float(match.group(6)))))
                    else:
                        #read cell
                        match = re.search('^.*VOLUME and BASIS-vectors are now.*$',line)
                        #" VOLUME and BASIS-vectors are now :"
                        if(match): # forces found
                            sublevel = "cell"
                            cell_frame   = []
                            continue
                        #read forces
                        match = re.search('^.*TOTAL-FORCE.*$',line)
                        # " POSITION                                       TOTAL-FORCE (eV/Angst)"
                        if(match): # forces found
                            sublevel = "forces"
                            positions_frame   = []
                            forces_frame   = []
                            continue
                        if(self.idipol > 0):
                            #read dipole moments
                            match = re.search('^[ \t]*dipolmoment[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*electrons.*$',line)
                            # " dipolmoment          -0.000002     -0.051875     -0.352631 electrons x Angstroem"
                            if(match): # dipol moment found
                                dipol_store = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                                if(self.ldipol == False):
                                    self.dipols.append(dipol_store)
                                continue
        
                            #for LDIPOL = T, dipols are found for every electronic step, so only use the converged step
                            if(self.ldipol == True):
                                match = re.search('^.*aborting loop because EDIFF is reached.*$',line)
                                # "------------------------ aborting loop because EDIFF is reached ----------------------------------------"
                                if(match): # electronic convergence
                                    self.dipol_diff.append(dipol_store)
                                    dipol_store = None
                                    continue
    
                        # Check for end of OUTCAR.
                        if(re.search('^[ \t]*General timing and accounting informations for this job:[ \t]*$',line)):
                            readmode += 1
                            continue
                elif(readmode == 4):
                    # print("self.forces_ref")
                    # for force in self.forces_ref:
                    #     print(force)
                    # print("self.forces_diff")
                    # for diff in self.forces_diff:
                    #     print("---")
                    #     for force in diff:
                    #         print(force)
                    break
        if(readmode == 0):
            print("Error: Could not find Vasp version")
        elif(readmode ==1):
            print("Error: Could not find elemental properties")
        elif(readmode ==2):
            print("Error: Could not find DIPOL tags")
        elif(readmode ==3):
            print("Error: Could not find end of OUTCAR")
    def getElement(self, number_in):
        atomsum = 0
        for number, element in zip(self.elnr, self.elements):
            atomsum += number
            if(number_in <= atomsum):
                return element
        print("Could not get element name")
        sys.exit()

    def getDistance(self, atom1, atom2, printmode = "all",pbc = False):
        if(printmode == "all"):
            for index in range(0,self.numframes):
                cart1 = self.positions[index][atom1-1]
                cart2 = self.positions[index][atom2-1]
                distances = np.dtype(float)
                distance_values = []
                distance_values.append(np.linalg.norm(cart2-cart1))
                if(pbc == True):
                    x = self.cell[index][0]
                    y = self.cell[index][1]
                    distance_values.append(np.linalg.norm(cart2-cart1+x))
                    distance_values.append(np.linalg.norm(cart2-cart1-x))
                    distance_values.append(np.linalg.norm(cart2-cart1+y))
                    distance_values.append(np.linalg.norm(cart2-cart1-y))
                    distance_values.append(np.linalg.norm(cart2-cart1+x+y))
                    distance_values.append(np.linalg.norm(cart2-cart1+x-y))
                    distance_values.append(np.linalg.norm(cart2-cart1-x+y))
                    distance_values.append(np.linalg.norm(cart2-cart1-x-y))
                distances = np.array(distance_values)
                print(distance.min())

        elif(printmode == "less"):
            cart1 = self.positions[-1][atom1-1]
            cart2 = self.positions[-1][atom2-1]
            distances = np.dtype(float)
            distance_values = []
            distance_values.append(np.linalg.norm(cart2-cart1))
            if(pbc == True):
                x = self.cell[-1][0]
                y = self.cell[-1][1]
                distance_values.append(np.linalg.norm(cart2-cart1+x))
                distance_values.append(np.linalg.norm(cart2-cart1-x))
                distance_values.append(np.linalg.norm(cart2-cart1+y))
                distance_values.append(np.linalg.norm(cart2-cart1-y))
                distance_values.append(np.linalg.norm(cart2-cart1+x+y))
                distance_values.append(np.linalg.norm(cart2-cart1+x-y))
                distance_values.append(np.linalg.norm(cart2-cart1-x+y))
                distance_values.append(np.linalg.norm(cart2-cart1-x-y))
            distances = np.array(distance_values)
            print(distances.min())
    def getDipols(self):
        pass
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
    def writeDipols(self,printmode = "all"):
        pass