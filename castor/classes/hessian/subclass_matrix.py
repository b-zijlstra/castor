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

#MY PATHS

#MY CLASSES
from subclass_frequency import Frequency

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
        if(self.nonsym == None or len(self.nonsym) == 0 or len(self.nonsym[0]) == 0):
            print "Can not symmetrize matrix because not symmetrized matrix does not exist!"
            sys.exit()
        else:
            self.sym = []
            for i in range(0,len(self.nonsym)):
                row = []
                for j in range(0,len(self.nonsym)):
                    if i == j:
                        row.append(self.nonsym[i][j])
                    else:
                        row.append(0.5*(self.nonsym[i][j]+self.nonsym[j][i]))
                self.sym.append(row)
    def sym2mass(self, map_in):
        if(self.sym == None or len(self.sym) == 0 or len(self.sym[0]) == 0):
            print "Can not make mass matrix because symmetrized matrix does not exist!"
            sys.exit()
        else:
            self.mass = []
            for row, atom_row in zip(self.sym, self.atoms):
                mass_row = map_in[atom_row]
                rowlist = []
                for column, atom_col in zip(row, self.atoms):
                    mass_col = map_in[atom_col]
                    rowlist.append(column / math.sqrt(mass_row*mass_col))
                self.mass.append(rowlist)
    def mass2diag(self):
        if(self.mass == None or len(self.mass) == 0 or len(self.mass[0]) == 0):
            print "Can not diagonalize mass matrix because it does not exist!"
            sys.exit()
        else:
            numpymat = np.matrix(self.mass)
            self.diag = np.linalg.eigh(numpymat)
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
            string += '{0:>{width}} '.format(count, width=self.spaces+1)
            string += freq.getString()
            print string
    def printMatrix(self, matrix_in):
        self.printAxis_X()
        for row, atom_row, axes_row in zip(matrix_in, self.atoms, self.axes):
            name = str(atom_row)+axes_row
            string = '{0:^{width}}'.format(name, width=self.spaces+2)
            for column in row:
                data = '{0:.{width}f}'.format(column, width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+3)
            print string
    def printAxis_X(self):
        string = '{0:^{width}}'.format("", width=self.spaces+2)
        for atom, axis in zip(self.atoms, self.axes):
            name = str(atom)+axis
            string += '{0:>{width}}'.format(name, width=self.decimals+self.spaces+3)
        print string