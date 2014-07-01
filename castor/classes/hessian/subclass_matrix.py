#! /usr/bin/env python

# 
# subclass_matrix.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math
import numpy as np
from copy import deepcopy as copy


#MY PATHS

#MY CLASSES
from classes.hessian.subclass_frequency import Frequency

class Matrix:
    """Defines the matrices of a Hessian"""
    def __init__(self, atoms_in = None, axes_in = None):
        self.decimals = 6
        self.spaces = 3
        self.atoms = atoms_in
        self.axes = axes_in
        self.frequencies = []
    def setup(self, nonsym_in = None, sym_in = None, mass_in = None, diag_in = None):
        self.nonsym = nonsym_in
        self.sym = sym_in
        self.mass = mass_in
        self.diag = diag_in
    def nonsym2sym(self):
        "nonsym2sym is not yet implemented"
        sys.exit()
    def sym2mass(self, map_in):
        self.mass = copy(self.sym)
        for row, atom_row in zip(mass, self.atoms):
            mass_row = self.massmap[atom_row]
            for column, atom_col in zip(row, self.atoms):
                mass_col = map_in[atom_col]
                column /= math.sqrt(mass_row*mass_col)
    def mass2diag(self):
        if(self.mass == None or len(self.mass) == 0 or self.mass[0].size == 0):
            print "Can not diagonalize mass matrix because it does not exist!"
            sys.exit()
        else:
            self.diag = np.linalg.eigh(self.mass)
    def diag2freq(self):
        for eigenval in self.diag[0]:
            if(eigenval < 0):
                eigenroot = math.sqrt(abs(eigenval))
                self.frequencies.append(Frequency(eigenroot,False))
            else:
                eigenroot = math.sqrt(abs(eigenval))
                self.frequencies.append(Frequency(eigenroot,True))
    def printNonsym(self):
        self.printMatrix(self.nonsym)
    def printSym(self):
        self.printMatrix(self.sym)
    def printMass(self):
        self.printMatrix(self.mass)
    def printDiag(self):
        self.printMatrix(self.diag)
    def printFreq(self):
        for freq, count in zip(self.frequencies, range(1,len(self.frequencies)+1)):
            string = ""
            string += '{:>{width}} '.format(count, width=self.spaces+1)
            string += freq.getString()
            print string
    def printMatrix(self, matrix_in):
        self.printAxis_X()
        for row, atom_row, axes_row in zip(matrix_in, self.atoms, self.axes):
            name = str(atom_row)+axes_row
            string = '{:^{width}}'.format(name, width=self.spaces+2)
            print row
            for column in row:
                print column
                data = '{:.{width}f}'.format(column, width=self.decimals)
                string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+3)
            print string
    def printAxis_X(self):
        string = '{:^{width}}'.format("", width=self.spaces+2)
        for atom, axis in zip(self.atoms, self.axes):
            name = str(atom)+axis
            string += '{:>{width}}'.format(name, width=self.decimals+self.spaces+3)
        print string