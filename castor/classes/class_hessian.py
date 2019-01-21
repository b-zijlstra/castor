#! /usr/bin/env python3

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
from .hessian.subclass_matrix import Matrix
from .hessian.subclass_matrix import Dipols
from .hessian.subclass_frequency import Frequency

class Hessian:
    """Defines a Hessian"""
    def __init__(self, temp_in = [800]):
        # General storage
        self.ibrion      = None             # IBRION tag (should be 5).
        self.decimals    = 6                # Number of decimals to output.
        self.spaces      = 3                # Number of spaces to output.
        self.elnr        = np.dtype(int)    # Vector of element numbers.
        self.elements    = []               # List of element names.
        self.masses      = []               # List of element masses.
        self.temp        = temp_in          # Temperature in Kelvin - Used for the vibrational partition function
        self.kb          = 8.617332478e-2   # Boltzmann constant in meV/K
        self.skipset     = set()            # Set of atoms that should be skipped.

        # Storage for Hessian matrix
        self.finitdiff   = None             # Finite difference used.
        self.forces_ref  = None             # Storage for the forces of the initial system.
        self.forces_diff = []               # List of forces for each displacement
        self.matrix      = None             # Hessian matrix container.
        self.numberset   = set()            # Set of atom numbers listed for modification.
        self.massmap     = dict()           # Map with masses corresponding to atoms.
        self.nummap      = dict()           # Map with modified masses corresponding to atoms.
        self.numbers     = None             # Optional selection of atoms.
        self.element     = None             # Selection of element to modify.
        self.mass        = None             # Mass to give to selected element.
        self.changes     = False            # Bool to log whether a change has been made.
        self.newmatrices = []               # Generated matrices with new masses.

        # Storage for Dipol matrix
        self.volume      = None             # Volume of the cell
        self.idipol      = None             # Which IDIPOL has been used.
        self.ldipol      = None             # Whether dipol corrections are included in the potential.
        self.epsilon     = None             # Bulk dielectric constant.
        self.dipols      = None             # Dipol matrix container.
        self.dipol_diff  = []               # List of dipol differences.
        self.frequencies = []               # List of dipol frequencies.
        self.skipfreq    = []               # Reduced list of dipol frequencies without skipped atoms.
        self.freedom     = []               # List of degrees of freedom.
        self.dipol_ref   = None             # Storage for the dipol moment of the initial system.



    def read(self, outcar_in, mode_in):
        if(mode_in == "calc"):
            recalc_hessian = True
            recalc_cor = True
        else:
            recalc_hessian = False
            recalc_cor = False
        readmode = 0
        sublevel = None
        with open(outcar_in, 'r') as inputfile:
            # OUTCAR info
            vaspversion = None   # Which Vasp version has been used.
            nwrite      = None   # Which nwrite has been used.
            self.idipol = None   # Which IDIPOL has been used.
            self.ldipol = None   # Whether dipol corrections are included in the potential.

            # Temporary storage for Hessian matrix information
            floatlist_masses           = []     # List of atom masses.
            floatlist_list_rows_nonsym = []     # List of nonsym matrix rows.
            floatlist_list_rows_sym    = []     # List of sym matrix rows.
            floatlist_list_rows_mass   = []     # List of mass matrix rows.
            scale_diff                 = False  # Whether displacements should be corrected for mass.
            
            # Temporary storage for Dipol matrix
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
                    # What element mass?
                    match = re.search('^[ \t]*POMASS[ \t]+=[ \t]+([0-9.]+);.*$',line)
                    if(match):
                        self.masses.append(float(match.group(1)))
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
                        if(self.ibrion != 5):
                            print("Error: IBRION != 5: " + line)
                            sys.exit()
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
                                if(self.idipol==0):
                                    recalc_cor = False
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
                                recalc_cor = False
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
                elif(readmode == 3): # Read the forces and read the dipolmoments if IDIPOL > 0
                    if(sublevel == "forces"):
                        match = re.search('^.*total drift:.*$',line)
                        # "    total drift:                               -0.002172     -0.000568      0.002532"
                        if(match): # all forces read
                            sublevel = None
                            if(self.forces_ref == None):
                                self.forces_ref = forces
                            else:
                                self.forces_diff.append(forces)
                            continue
                        match = re.search('^[ \t]*([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]*$',line)
                        # "      0.72358      1.25327      9.49588         0.001132      0.001590      0.000119"
                        if(match): # force found
                            forces.append([float(match.group(4)), float(match.group(5)), float(match.group(6))])
                    else:
                        #read forces
                        match = re.search('^.*TOTAL-FORCE.*$',line)
                        # " POSITION                                       TOTAL-FORCE (eV/Angst)"
                        if(match): # forces found
                            sublevel = "forces"
                            forces   = []
                            continue
                        if(self.idipol > 0):
                            #read dipole moments
                            match = re.search('^[ \t]*dipolmoment[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*([0-9.-]+)[ \t]*electrons.*$',line)
                            # " dipolmoment          -0.000002     -0.051875     -0.352631 electrons x Angstroem"
                            if(match): # dipol moment found
                                dipol_store = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                                if(self.ldipol == False):
                                    if(self.dipol_ref==None):
                                        self.dipol_ref = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                                    else:
                                        dipol = [float(match.group(1)), float(match.group(2)), float(match.group(3))]
                                        self.dipol_diff.append(dipol)
                                continue
        
                            #for LDIPOL = T, dipols are found for every electronic step, so only use the converged step
                            if(self.ldipol == True):
                                match = re.search('^.*aborting loop because EDIFF is reached.*$',line)
                                # "------------------------ aborting loop because EDIFF is reached ----------------------------------------"
                                if(match): # electronic convergence
                                    if(self.dipol_ref==None):
                                        self.dipol_ref = dipol_store
                                        dipol_store    = None
                                    else:
                                        self.dipol_diff.append(dipol_store)
                                        dipol_store = None
                                    continue
    
                        # Check for not symmetrized hessian matrix.
                        if(re.search('^[ \t]*SECOND DERIVATIVES \(NOT SYMMETRIZED\)[ \t]*$',line)):
                            readmode += 1
                            continue
                elif(readmode == 4): # Read the degrees of freedom
                    # "              11X         11Y         11Z         12X         12Y         12Z"
                    match = re.search('^([ \t]+[0-9]+[XYZ])+[ \t]*$',line)
                    if(match):
                        self.freedom = match.group(0).split()
                        freedom_count = 0
                        readmode      += 1
                        continue
                elif(readmode == 5): # Read not symmetrized hessian matrix.
                    if(recalc_hessian == False):
                        # Filling the matrix
                        match = re.search('^[ \t]*([0-9]+)([XYZ])[ \t]*(.*[0-9])[ \t]*$',line)
                        if(match):
                            atom = int(self.freedom[freedom_count][:-1])
                            direction = self.freedom[freedom_count][-1]
                            if(atom != int(match.group(1)) or direction != match.group(2)):
                                print("Error: Degrees of freedom do not match: " + self.freedom[freedom_count] + " != " + match.group(1) + match.group(2))
                                sys.exit()
                            stringlist_row = match.group(3).split()
                            floatlist_values = []
                            for i in range(0,len(stringlist_row)):
                                floatlist_values.append(float(stringlist_row[i]))
                            floatlist_list_rows_nonsym.append(floatlist_values)
                            freedom_count += 1
                        if(freedom_count == len(self.freedom)):
                            matrix_size = len(floatlist_list_rows_nonsym)
                            self.matrix = Matrix(self.freedom)
                            # self.matrix.setup(rows_nonsym, rows_sym, rows_mass, diag)
                            self.matrix.setup(floatlist_list_rows_nonsym, None, None, None)
                            self.matrix.nonsym2sym()
                            self.mapMass()
                            # self.skipset is always an empty set here
                            self.matrix.sym2mass(self.massmap,self.skipset)
                            readmode += 1
                            continue
                    elif(recalc_hessian == True):
                        self.matrix = Matrix(self.freedom)
                        self.matrix.setup(None, None, None, None)
                        if(recalc_cor == True):
                            self.matrix.getnonsym(self.freedom, self.finitdiff, self.forces_diff, self.epsilon, self.volume, self.dipol_diff)
                        else:
                            self.matrix.getnonsym(self.freedom, self.finitdiff, self.forces_diff, self.epsilon, self.volume, [])
                        self.matrix.nonsym2sym()
                        self.mapMass()
                        # self.skipset is always an empty set here
                        self.matrix.sym2mass(self.massmap,self.skipset)
                        readmode += 1
                        continue
                elif(readmode == 6 or readmode == 7): # Search for dynamical matrix.
                    # Make sure we are reading the correct dynamical matrix
                    if(re.search('^[ \t]*Eigenvectors and eigenvalues of the dynamical matrix[ \t]*$',line)):
                        readmode += 1
                        if(vaspversion == 5 and nwrite < 3):
                            scale_diff = True
                            readmode   += 1
                        continue
                elif(readmode == 8):
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
                            self.frequencies[-1].atdiff[atnum] = diff
                            continue
                        match = re.search('^[ \t]*Finite differences POTIM=[ \t]*[0-9.E-]+[ \t]*$',line)
                        # "Finite differences POTIM=  2.00000000000000004E-002"
                        if(match): # ending frequency reading
                            readmode += 1 # making sure no more freqs are read
                            continue
                elif(readmode == 9):
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
            print("Error: Could not find start of hessian")
        elif(readmode ==4):
            print("Error: Could not read degrees of freedom")
        elif(readmode ==5):
            print("Error: Problem reading not symmetrized hessian matrix")
        elif(readmode ==6 or readmode == 7):
            print("Error: Could not find dynamical matrix")
        elif(readmode ==8):
            print("Error: Could not find all frequencies")
    
    def addmatrix(self):
        new = Matrix(self.matrix.freedom)
        new.setup(self.matrix.nonsym, self.matrix.sym, None, None)
        self.newmatrices.append(new)
    def mapMass(self, element_in = None, mass_in = None, numbers_in = None):
        self.element   = element_in
        self.mass      = mass_in
        self.numbers   = numbers_in
        self.numberset = self.string2numberset(numbers_in, element_in)
        
        self.massmap.clear()
        for degree in self.matrix.freedom:
            atom = int(degree[:-1])
            if(atom not in self.massmap):
                if(atom in self.numberset and mass_in != None):
                    if(self.changes == False and mass_in != self.getMass(atom)):
                        self.changes = True
                    self.massmap[atom] = mass_in
                else:
                    self.massmap[atom] = self.getMass(atom)
    def mapMass_nummap(self, nummap_in = None):
        self.nummap = self.string2numbermap(nummap_in)
        if(self.nummap == None):
            print("Error: nummap == None")
            sys.exit()
        
        self.massmap.clear()
        for degree in self.matrix.freedom:
            atom = int(degree[:-1])
            if(atom not in self.massmap):
                if(atom in self.nummap):
                    if(self.changes == False and self.nummap[atom] != self.getMass(atom)):
                        self.changes = True
                    self.massmap[atom] = self.nummap[atom]
                else:
                    self.massmap[atom] = self.getMass(atom)
    def setSkip(self, numbers_in = None):
        self.skipset = self.string2numberset(numbers_in, "all")
        if (len(self.skipset) > 0):
            self.changes = True
    def string2numberset(self, string_in, element_in):
        numberset = set()
        if(element_in != None):
            if(string_in == None or string_in == "all"):
                for degree in self.matrix.freedom:
                    atom = int(degree[:-1])
                    if(self.getElement(atom)==element_in):
                        numberset.add(atom)
                return numberset
            else:
                maxatom = 0
                for degree in self.matrix.freedom:
                    atom = int(degree[:-1])
                    if(atom>maxatom):
                        maxatom = atom
                templist = string_in.split(',')
                templist = [x.strip() for x in templist]
                skiporallow = "skip"
                for x in templist:
                    if(x == "a"):
                        skiporallow = "allow"
                    else:
                        try:
                            x = [int(x)]
                        except ValueError:
                            x = x.split('-')
                            try:
                                x = range(int(x[0]),int(x[1])+1)
                            except ValueError:
                                try:
                                    if(x[1].strip()==":"):
                                        x = range(int(x[0]),maxatom+1)
                                    else:
                                        print("Could not set numbers setting")
                                        sys.exit()
                                except ValueError:
                                    print("Could not set numbers setting")
                                    sys.exit()
                        for y in x:
                            if(element_in == "all" or self.getElement(y)==element_in):
                                numberset.add(y)
                if(skiporallow == "skip"):
                    return numberset
                else:
                    skipset = set()
                    for degree in self.matrix.freedom:
                        atom = int(degree[:-1])
                        if(atom in numberset):
                            pass
                        else:
                            skipset.add(atom)
                    return skipset
        else:
            return numberset
    def string2numbermap(self, string_in):
        numbermap = dict()
        maxatom = 0
        for degree in self.matrix.freedom:
            atom = int(degree[:-1])
            if(atom>maxatom):
                maxatom = atom
        templist = string_in.split(',')
        templist = [x.strip() for x in templist]
        for x in templist:
            pair = x.split('=')
            try:
                pair[0] = int(pair[0])
                pair[1] = float(pair[1])
            except ValueError:
                print("Could not set numbers setting for pair:")
                print(pair)
                sys.exit()
            numbermap[pair[0]] = pair[1]
        return numbermap
    def getElement(self, number_in):
        atomsum = 0
        for number, element in zip(self.elnr, self.elements):
            atomsum += number
            if(number_in <= atomsum):
                return element
        print("Could not get element name")
        sys.exit()
    def getMass(self, number_in):
        atomsum = 0
        for number, mass in zip(self.elnr, self.masses):
            atomsum += number
            if(number_in <= atomsum):
                return mass
        print("Could not get element mass")
        sys.exit()
    def calcZPE(self, freqs_in):
        zpe = 0.0
        for freq in freqs_in:
            if(freq.imaginary==False):
                zpe += freq.meV
        return 0.0005 * zpe
    def calcPartition(self, freqs_in, temp_in):
        kbT = temp_in * self.kb
        nu  = 1.0
        for freq in freqs_in:
            if(freq.imaginary==False):
                nu *= 1.0 / (1.0 - math.exp(-freq.meV / kbT))
        return nu
    def write(self, printmode,first=None):
        if(printmode == "normal" or printmode == "all"):
            print("---------------------------")
            print("-        Settings:        -")
            print("---------------------------")
            if(len(self.nummap)>0):
                print("Using input nummap to set new mass.")
            elif(self.element == None):
                print("No elements selected to set new mass.")
            else:
                print("Element to set new mass = " + str(self.element))
                if(self.numbers == None):
                    print("Atom numbers of element to set new mass = all")
                else:
                    print("Atom numbers of element to set new mass = " + str(self.numbers))
                print("New mass = " + str(self.mass))
            print("\n")
            if(self.idipol > 0):
                print("---------------------------")
                print("-   Dipole information:   -")
                print("---------------------------")
                print("IDIPOL = " + str(self.idipol))
                print("Initial dipole moment: " + self.writeList3(self.dipol_ref))
                self.dipols.printMatrix()
                print("\n")
            if(printmode == "all" and self.matrix != None):
                # print("---------------------------")
                # print("-  Original nonsym matrix:  -")
                # print("---------------------------")
                # self.matrix.printNonsym()
                # print("---------------------------")
                # print("-  Original sym matrix:  -")
                # print("---------------------------")
                # self.matrix.printSym()
                print("---------------------------")
                print("-  Original mass matrix:  -")
                print("---------------------------")
                self.matrix.printMass()
                print("\n")
            for newmatrix in self.newmatrices:
                print("---------------------------")
                print("-        Changes:         -")
                print("---------------------------")
                self.printChanges()
                print("\n")
                if(printmode == "all"):
                    print("---------------------------")
                    print("-    New mass matrix:     -")
                    print("---------------------------")
                    newmatrix.printMass()
                    print("\n")
            if(self.matrix != None):
                print("---------------------------")
                print("-  Original frequencies:  -")
                print("---------------------------")
                # self.matrix.printFreq()
                count = 0
                for freq in self.matrix.frequencies:
                    count += 1
                    if(first != None and count > first):
                        break
                    freq.write(prefix=count, printmode=printmode)
                    if(printmode=="all"):
                        freq.writeDiff()
                print("\n")
                zpe = self.calcZPE(self.matrix.frequencies)
                print('Total ZPE contribution of frequencies in eV: {0:.{width}f}'.format(zpe, width=self.decimals))
                for temperature in self.temp:
                    nu = self.calcPartition(self.matrix.frequencies,temperature)
                    print('Vibrational partition function at {0:.2f} K: {1:.{width}E}'.format(temperature,nu, width=self.decimals))
                print("\n")
            for newmatrix in self.newmatrices:
                print("---------------------------")
                print("-    New frequencies:     -")
                print("---------------------------")
                # newmatrix.printFreq()
                count = 0
                for freq in newmatrix.frequencies:
                    count += 1
                    if(first != None and count > first):
                        break
                    freq.write(prefix=count, printmode=printmode)
                    if(printmode=="all"):
                        freq.writeDiff()
                print("\n")
                zpe = self.calcZPE(newmatrix.frequencies)
                print('Total ZPE contribution of frequencies in eV: {0:.{width}f}'.format(zpe, width=self.decimals))
                for temperature in self.temp:
                    nu = self.calcPartition(self.newmatrix.frequencies,temperature)
                    print('Vibrational partition function at {0:.2f} K: {1:.{width}E}'.format(temperature,nu, width=self.decimals))
                print("\n")
        elif(self.matrix != None):
            zpe = self.calcZPE(self.matrix.frequencies)
            if(printmode == "least"):
                print('{0:.{width}f}'.format(zpe, width=self.decimals))
                for temperature in self.temp:
                    nu = self.calcPartition(self.matrix.frequencies,temperature)
                    print('{0:.{width}E}'.format(nu, width=self.decimals))
            else:
                print('Total ZPE contribution of frequencies in eV: {0:.{width}f}'.format(zpe, width=self.decimals))
                for temperature in self.temp:
                    nu = self.calcPartition(self.matrix.frequencies,temperature)
                    print('Vibrational partition function at {0:.2f} K: {1:.{width}E}'.format(temperature,nu, width=self.decimals))
                if(self.idipol > 0):
                    count    = 0
                    highfreq = None
                    highint  = None
                    for freq in self.matrix.frequencies:
                        count += 1
                        if(count == 1):
                            highfreq = freq
                            highint  = freq
                            freqcount = count
                            intcount = count
                        else:
                            if(freq.cm1 > highint.cm1):
                                highfreq = freq
                                freqcount = count
                            if(freq.intensity[3] > highint.intensity[3]):
                                highint = freq
                                intcount = count
                    hf_f = highfreq.cm1
                    hf_i = highfreq.intensity[3]
                    hi_f = highint.cm1
                    hi_i = highint.intensity[3]
                    hf_str = '{0:>{width}}'.format(freqcount, width=self.spaces)
                    hf_str += " "
                    if(highfreq.imaginary==False):
                        hf_str += "f  ="
                    elif(highfreq.imaginary==True):
                        hf_str += "f/i="
                    hi_str = '{0:>{width}}'.format(intcount, width=self.spaces)
                    hi_str += " "
                    if(highint.imaginary==False):
                        hi_str += "f  ="
                    elif(highint.imaginary==True):
                        hi_str += "f/i="
                    print('Highest frequency: {0} {1:.{width}f} cm-1 | {2:.{width}f} int'.format(hf_str, hf_f,hf_i, width=self.decimals))
                    print('Highest intensity: {0} {1:.{width}f} cm-1 | {2:.{width}f} int'.format(hi_str, hi_f,hi_i, width=self.decimals))
            for newmatrix in self.newmatrices:
                zpe = self.calcZPE(newmatrix.frequencies)
                if(printmode == "least"):
                    print('{0:.{width}f}'.format(zpe, width=self.decimals))
                    for temperature in self.temp:
                        nu = self.calcPartition(newmatrix.frequencies,temperature)
                        print('{0:.{width}E}'.format(nu, width=self.decimals))
                else:
                    self.printChanges()
                    print('Total ZPE contribution of new frequencies in eV: {0:.{width}f}'.format(zpe, width=self.decimals))
                    for temperature in self.temp:
                        nu = self.calcPartition(newmatrix.frequencies,temperature)
                        print('New vibrational partition function at {0:.2f} K: {1:.{width}E}'.format(temperature,nu, width=self.decimals))
                    if(self.idipol > 0):
                        count    = 0
                        highfreq = None
                        highint  = None
                        for freq in newmatrix.frequencies:
                            count += 1
                            if(count == 1):
                                highfreq = freq
                                highint  = freq
                                freqcount = count
                                intcount = count
                            else:
                                if(freq.cm1 > highint.cm1):
                                    highfreq = freq
                                    freqcount = count
                                if(freq.intensity[3] > highint.intensity[3]):
                                    highint = freq
                                    intcount = count
                        hf_f = highfreq.cm1
                        hf_i = highfreq.intensity[3]
                        hi_f = highint.cm1
                        hi_i = highint.intensity[3]
                        hf_str = '{0:>{width}}'.format(freqcount, width=self.spaces)
                        hf_str += " "
                        if(highfreq.imaginary==False):
                            hf_str += "f  ="
                        elif(highfreq.imaginary==True):
                            hf_str += "f/i="
                        hi_str = '{0:>{width}}'.format(intcount, width=self.spaces)
                        hi_str += " "
                        if(highint.imaginary==False):
                            hi_str += "f  ="
                        elif(highint.imaginary==True):
                            hi_str += "f/i="
                        print('Highest frequency: {0} {1:.{width}f} cm-1 | {2:.{width}f} int'.format(hf_str, hf_f,hf_i, width=self.decimals))
                        print('Highest intensity: {0} {1:.{width}f} cm-1 | {2:.{width}f} int'.format(hi_str, hi_f,hi_i, width=self.decimals))
    def printChanges(self):
        if(len(self.numberset)>0):
            for i in self.numberset:
                string = "El: "
                string += self.getElement(i)
                string += "  Nr: "
                string += str(i)
                string += "    Mass: "
                string += str(self.getMass(i))
                string += " -> "
                string += str(self.massmap[i])
                print(string)
        if(len(self.skipset)>0):
            for i in self.skipset:
                string = "Skip atom nr: "
                string += str(i)
                print(string)
        if(len(self.nummap)>0):
            for i in self.nummap:
                string = "El: "
                string += self.getElement(i)
                string += "  Nr: "
                string += str(i)
                string += "    Mass: "
                string += str(self.getMass(i))
                string += " -> "
                string += str(self.massmap[i])
                print(string)
        if(len(self.numberset)==0 and len(self.skipset)==0 and len(self.nummap)==0):
            print("None")
    def getDipols(self, frequencies_in = None):
        self.dipols = Dipols()
        self.dipols.setup(self.dipol_diff,self.freedom,self.finitdiff)
        self.dipols.setIntensities(self.frequencies)
        self.getIntens(self.matrix.frequencies)
    def getIntens(self,frequencies_in):
        self.dipols.setIntensities(frequencies_in)
    def getNewIntens(self):
        for newmat in self.newmatrices:
            self.getIntens(newmat.frequencies)
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
        if(printmode == "normal" or printmode=="all"):
            print("IDIPOL = " + str(self.idipol))
            print("Initial dipole moment: " + self.writeList3(self.dipol_ref))
            self.dipols.printMatrix()
        count = 0
        for freq in self.frequencies:
            count += 1
            freq.write(prefix=count, printmode=printmode)