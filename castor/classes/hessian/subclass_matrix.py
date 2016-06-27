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
    def __init__(self, freedom_in):
        self.decimals    = 6
        self.spaces      = 3
        self.freedom     = freedom_in
        self.frequencies = []
        self.skipset     = set()
    def setup(self, nonsym_in = None, sym_in = None, mass_in = None, diag_in = None):
        self.nonsym = nonsym_in
        self.sym    = sym_in
        self.mass   = mass_in
        self.diag   = diag_in
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
            for row, degree in zip(self.sym, self.freedom):
                atom     = int(degree[:-1])
                mass_row = map_in[atom]
                rowlist  = []
                if(atom in self.skipset):
                    redefineAxes = True # Unused
                    continue
                for column, degree in zip(row, self.freedom):
                    atom     = int(degree[:-1])
                    mass_col = map_in[atom]
                    if(atom in self.skipset):
                        redefineAxes = True # Unused
                    else:
                        rowlist.append(column / math.sqrt(mass_row*mass_col))
                self.mass.append(rowlist)
    def mass2diag(self):
        # self.printMass()
        if(self.mass == None or len(self.mass) == 0 or len(self.mass[0]) == 0):
            print "Can not diagonalize mass matrix because it does not exist!"
            sys.exit()
        else:
            numpymat = np.matrix(self.mass)
            self.diag = np.linalg.eigh(numpymat)
    def diag2freq(self, map_in):
        freedom_reduced     = []
        set_reduced = set()
        for degree in self.freedom:
            atom = int(degree[:-1])
            if(atom not in self.skipset):
                set_reduced.add(atom)
                freedom_reduced.append(degree)

        # fulldiff will contain dict() for every DOF/FREQ
        # dict() will contain atom number and [x.x, x.x, x.x] diff
        fulldiff=[]
        for degree in freedom_reduced:
            dictionary = dict()
            for atom in set_reduced:
                dictionary[atom] = [0.0,0.0,0.0]
            fulldiff.append(dictionary)

        for degree,diffvec in zip(freedom_reduced,self.diag[1].tolist()):
            for atdiff, diff in zip(fulldiff,diffvec):
                atom      = int(degree[:-1])
                direction = degree[-1]
                mass      = map_in[atom]
                if(direction == "X"):
                    atdiff[atom][0] = diff/math.sqrt(mass)
                elif(direction == "Y"):
                    atdiff[atom][1] = diff/math.sqrt(mass)
                elif(direction == "Z"):
                    atdiff[atom][2] = diff/math.sqrt(mass)
        # print "Filled fulldiff"
        # print fulldiff

        for eigenval,atdiff in zip(self.diag[0],fulldiff):
            newfreq = Frequency()
            eigenroot = math.sqrt(abs(eigenval))
            if(eigenval < 0):
                newfreq.setEigen(eigenroot,atdiff,False)
            elif(eigenval > 0):
                newfreq.setEigen(eigenroot,atdiff,True)
            self.frequencies.append(newfreq)
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
        if(len(matrix_in) < len(self.freedom)):
            self.printAxis_X(True)
            freedomlist = []
            for degree in self.freedom:
                atom = int(degree[:-1])
                if (atom in self.skipset):
                    pass
                else:
                    freedomlist.append(degree)
        else:
            self.printAxis_X(False)
            freedomlist = self.freedom
        for row, degree in zip(matrix_in, freedomlist):
            string = '{0:^{width}}'.format(degree, width=self.spaces+2)
            for column in row:
                data = '{0:.{width}f}'.format(column, width=self.decimals)
                string += '{0:>{width}}'.format(data, width=self.decimals+self.spaces+3)
            print string
    def printAxis_X(self,skip_bool = False):
        string = '{0:^{width}}'.format("", width=self.spaces+2)
        for degree in self.freedom:
            atom = int(degree[:-1])
            if (skip_bool == True and atom in self.skipset):
                pass
            else:
                string += '{0:>{width}}'.format(degree, width=self.decimals+self.spaces+3)
        print string

class Dipols:
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
            sumx = 0
            sumy = 0
            sumz = 0
            for dictionary in freq.atdiff:
                atom = dictionary
                diff = freq.atdiff[dictionary]
                for dip in self.dipols:
                    if(dip[0] == atom):
                        if(dip[1] == "X"):
                            sumx+=dip[2]* diff[0]
                            sumy+=dip[3]* diff[0]
                            sumz+=dip[4]* diff[0]
                        elif(dip[1] == "Y"):
                            sumx+=dip[2]* diff[1]
                            sumy+=dip[3]* diff[1]
                            sumz+=dip[4]* diff[1]
                        elif(dip[1] == "Z"):
                            sumx+=dip[2]* diff[2]
                            sumy+=dip[3]* diff[2]
                            sumz+=dip[4]* diff[2]
            totaldipol     = sumx**2 + sumy**2 + sumz**2
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