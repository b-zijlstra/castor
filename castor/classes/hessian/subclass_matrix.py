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
from subclass_frequency import DipolFrequency

class Matrix:
    """Defines the matrices of a Hessian"""
    def __init__(self, atoms_in = None, axes_in = None):
        self.decimals = 6
        self.spaces = 3
        self.atoms = atoms_in
        self.axes = axes_in
        self.frequencies = []
        self.skipset = set()
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
    def sym2mass(self, map_in,set_in):
        self.skipset = set_in
        if(self.sym == None or len(self.sym) == 0 or len(self.sym[0]) == 0):
            print "Can not make mass matrix because symmetrized matrix does not exist!"
            sys.exit()
        else:
            self.mass = []
            for row, atom_row in zip(self.sym, self.atoms):
                mass_row = map_in[atom_row]
                rowlist = []
                if(atom_row in set_in):
                    redefineAxes = True
                    continue
                for column, atom_col in zip(row, self.atoms):
                    mass_col = map_in[atom_col]
                    if(atom_col in set_in):
                        redefineAxes = True
                    else:
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
            elif(eigenval > 0):
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
        if(len(matrix_in) < len(self.atoms)):
            self.printAxis_X(True)
            atomlist = []
            axeslist = []
            for atom, axis in zip(self.atoms, self.axes):
                if (atom in self.skipset):
                    pass
                else:
                    atomlist.append(atom)
                    axeslist.append(axis)
        else:
            self.printAxis_X(False)
            atomlist = self.atoms
            axeslist = self.axes
        for row, atom_row, axes_row in zip(matrix_in, atomlist, axeslist):
            name = str(atom_row)+axes_row
            string = '{0:^{width}}'.format(name, width=self.spaces+2)
            for column in row:
                data = '{0:.{width}f}'.format(column, width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+3)
            print string
    def printAxis_X(self,skip_bool = False):
        string = '{0:^{width}}'.format("", width=self.spaces+2)
        for atom, axis in zip(self.atoms, self.axes):
            if (skip_bool == True and atom in self.skipset):
                pass
            else:
                name = str(atom)+axis
                string += '{0:>{width}}'.format(name, width=self.decimals+self.spaces+3)
        print string

class Dipols:
    """Defines the dipol matrix from a Hessian"""
    def __init__(self):
        self.decimals = 6
        self.spaces = 3
        self.dipols = []
    def setup(self, dipol_diff, degrees_freedom, diff):
        #check if the amount of dipols match the degrees of freedom
        if(len(dipol_diff)!= 2*len(degrees_freedom)):
            print "Error, amounts of dipols and degrees of freedom do not match!"
            sys.exit()

        atom = 0
        direction = 0
        dipcount = 0
        for degree in degrees_freedom:
            atom = int(degree[:-1])
            direction = degree[-1]
            #first dipol_diff is +diff, second dipol_diff is -diff
            mux = (dipol_diff[dipcount][0]-dipol_diff[dipcount+1][0])/(2*diff)
            muy = (dipol_diff[dipcount][1]-dipol_diff[dipcount+1][1])/(2*diff)
            muz = (dipol_diff[dipcount][2]-dipol_diff[dipcount+1][2])/(2*diff)
            self.dipols.append([atom,direction,mux, muy, muz])
            dipcount += 2
    def setIntensities(self,frequencies):
        for freq in frequencies:
            sumx=0
            sumy=0
            sumz=0
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
            xcor = '{0:.{width}f}'.format(dipol[2], width=self.decimals)
            ycor = '{0:.{width}f}'.format(dipol[3], width=self.decimals)
            zcor = '{0:.{width}f}'.format(dipol[4], width=self.decimals)
            string += '{0:>{width}}'.format(xcor, width=self.decimals+self.spaces+4)
            string += '{0:>{width}}'.format(ycor, width=self.decimals+self.spaces+4)
            string += '{0:>{width}}'.format(zcor, width=self.decimals+self.spaces+4)
            print string