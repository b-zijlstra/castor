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
from hessian.subclass_matrix import Matrix

class Hessian:
    """Defines a Hessian"""
    def __init__(self, temp_in = 800):
        self.decimals = 6
        self.spaces = 3
        self.matrix = None
        self.elnr = np.dtype(int)
        self.elements = np.dtype(str)
        self.masses = np.dtype(float)
        self.temp = temp_in # in Kelvin - Used for the vibrational partition function
        self.kb = 8.617332478e-2 # Botzmann constant in meV/K
        self.numberset = set()
        self.massmap = dict()
        self.numbers = None
        self.element = None
        self.mass = None
        self.changes = False
        self.newmatrices = []
    def read(self, outcar_in):
        with open(outcar_in, 'r') as inputfile:
            intlist_atoms = []
            charlist_axes = []
            stringlist_elements = []
            floatlist_masses = []
            floatlist_list_rows = []
            lines = inputfile.readlines()
            for line in lines:
                match = re.search('^[ \t]*([0-9]+)([XYZ])[ \t]*(.*[0-9])[ \t]*$',line);
                if(match):
                    intlist_atoms.append(int(match.group(1)))
                    charlist_axes.append(match.group(2))
                    stringlist_row = match.group(3).split()
                    floatlist_values = []
                    for i in range(0,len(stringlist_row)):
                        floatlist_values.append(float(stringlist_row[i]))
                    floatlist_list_rows.append(floatlist_values)
                match = re.search('^[ \t]*TITEL[ \t]+=[ \t]+[A-Z_]+ ([A-Za-z]+) [0-9A-Za-z]+[ \t]*$',line);
                if(match):
                    stringlist_elements.append(match.group(1))
                match = re.search('^[ \t]*POMASS[ \t]+=[ \t]+([0-9.]+);.*$',line);
                if(match):
                    floatlist_masses.append(float(match.group(1)))
                match = re.search('^[ \t]*ions per type =([0-9 \t]+)[ \t]*$',line);
                if(match):
                    stringlist_values = match.group(1).split()
                    intlist_values = []
                    for i in range(0,len(stringlist_values)):
                        intlist_values.append(int(stringlist_values[i]))
                    self.vector_elnr = np.array(intlist_values)
            matrix_size = len(floatlist_list_rows)/3
            self.matrix = Matrix(intlist_atoms[0:matrix_size], charlist_axes[0:matrix_size])
            self.matrix.setup(floatlist_list_rows[0 : matrix_size], floatlist_list_rows[matrix_size : matrix_size*2], floatlist_list_rows[matrix_size*2 : matrix_size*3], None)
            self.vector_elements = stringlist_elements
            self.vector_masses = floatlist_masses
        if(self.matrix.sym == None or len(self.matrix.sym) == 0 or len(self.matrix.sym[0]) == 0):
            print "Reading from " + outcar_in + " failed! Could not find symmetry matrix."
            sys.exit()
        if(self.matrix.mass == None or len(self.matrix.mass) == 0 or len(self.matrix.mass[0]) == 0):
            print "Reading from " + outcar_in + " failed! Could not find mass matrix."
            sys.exit()
        if(self.vector_elnr == None):
            print "Reading from " + outcar_in + " failed! Could not find number of elements."
            sys.exit()
    def addmatrix(self):
        new = Matrix(self.matrix.atoms, self.matrix.axes)
        new.setup(self.matrix.nonsym, self.matrix.sym, None, None)
        self.newmatrices.append(new)
    def mapMass(self, element_in = None, mass_in = None, numbers_in = None):
        self.element = element_in
        self.mass = mass_in
        self.numbers = numbers_in

        self.numberset = self.string2numberset(numbers_in, element_in)
        
        self.massmap.clear()
        for atom in self.matrix.atoms:
            if(atom not in self.massmap):
                if(atom in self.numberset and mass_in != None):
                    if(self.changes == False and mass_in != self.getMass(atom)):
                        self.changes = True
                    self.massmap[atom] = mass_in
                else:
                    self.massmap[atom] = self.getMass(atom)
    def string2numberset(self, string_in, element_in):
        numberset = set()
        if(string_in == None or string_in == "all"):
            for atom in self.matrix.atoms:
                if(self.getElement(atom)==element_in):
                    numberset.add(atom)
            return numberset
        else:
            templist = string_in.split(',')
            templist = [x.strip() for x in templist]
            for x in templist:
                try:
                    x = [int(x)]
                except ValueError:
                    x = x.split('-')
                    try:
                        x = range(int(x[0]),int(x[1])+1)
                    except ValueError:
                        try:
                            if(x[1].strip()==":"):
                                x = range(int(x[0]),np.max(self.matrix.atoms)+1)
                            else:
                                print "Could not set numbers setting"
                                sys.exit()
                        except ValueError:
                            print "Could not set numbers setting"
                            sys.exit()
                for y in x:
                    if(self.getElement(y)==element_in):
                        numberset.add(y)
            return numberset
    def getElement(self, number_in):
        atomsum = 0
        for number, element in zip(self.vector_elnr, self.vector_elements):
            atomsum += number
            if(number_in <= atomsum):
                return element
        print "Could not get element name"
        sys.exit()
    def getMass(self, number_in):
        atomsum = 0
        for number, mass in zip(self.vector_elnr, self.vector_masses):
            atomsum += number
            if(number_in <= atomsum):
                return mass
        print "Could not get element mass"
        sys.exit()
    def calcZPE(self, freqs_in):
        zpe = 0.0
        for freq in freqs_in:
            if(freq.imaginary==False):
                zpe += freq.meV
        return 0.0005 * zpe
    def calcPartition(self, freqs_in):
        kbT = self.temp * self.kb
        nu = 1.0
        for freq in freqs_in:
            if(freq.imaginary==False):
                nu *= 1.0 / (1.0 - math.exp(-freq.meV / kbT))
        return nu
    def write(self, printmode):
        if(printmode == "all"):
            print "---------------------------"
            print "-        Settings:        -"
            print "---------------------------"
            print "Element to set new mass = " + str(self.element)
            if(self.numbers == None):
                print "Atom numbers of element to set new mass = all"
            else:
                print "Atom numbers of element to set new mass = " + str(self.numbers)
            print "New mass = " + str(self.mass)
            print "\n"
            if(self.matrix != None):
                print "---------------------------"
                print "-  Original mass matrix:  -"
                print "---------------------------"
                self.matrix.printMass()
                print "\n"
            for newmatrix in self.newmatrices:
                print "---------------------------"
                print "-        Changes:         -"
                print "---------------------------"
                self.printChanges()
                print "\n"
                print "---------------------------"
                print "-    New mass matrix:     -"
                print "---------------------------"
                newmatrix.printMass()
                print "\n"
            if(self.matrix != None):
                print "---------------------------"
                print "-  Original frequencies:  -"
                print "---------------------------"
                self.matrix.printFreq()
                print "\n"
                zpe = self.calcZPE(self.matrix.frequencies)
                nu = self.calcPartition(self.matrix.frequencies)
                print 'Total ZPE contribution of frequencies in eV: {0:.{width}f}'.format(zpe, width=self.decimals)
                print 'Vibrational partition function: {0:.{width}f}'.format(nu, width=self.decimals)
                print "\n"
            for newmatrix in self.newmatrices:
                print "---------------------------"
                print "-    New frequencies:     -"
                print "---------------------------"
                newmatrix.printFreq()
                print "\n"
                zpe = self.calcZPE(newmatrix.frequencies)
                nu = self.calcPartition(newmatrix.frequencies)
                print 'Total ZPE contribution of frequencies in eV: {0:.{width}f}'.format(zpe, width=self.decimals)
                print 'Vibrational partition function: {0:.{width}f}'.format(nu, width=self.decimals)
                print "\n"
        elif(self.matrix != None):
            zpe = self.calcZPE(self.matrix.frequencies)
            nu = self.calcPartition(self.matrix.frequencies)
            print 'Total ZPE contribution of frequencies in eV: {0:.{width}f}'.format(zpe, width=self.decimals)
            print 'Vibrational partition function: {0:.{width}f}'.format(nu, width=self.decimals)

            for newmatrix in self.newmatrices:
                self.printChanges()
                zpe = self.calcZPE(newmatrix.frequencies)
                nu = self.calcPartition(newmatrix.frequencies)
                print 'Total ZPE contribution of new frequencies in eV: {0:.{width}f}'.format(zpe, width=self.decimals)
                print 'New vibrational partition function: {0:.{width}f}'.format(nu, width=self.decimals)
    def printChanges(self):
        if(len(self.numberset)==0):
            print "None"
        else:
            for i in self.numberset:
                string = "El: "
                string += self.getElement(i)
                string += "  Nr: "
                string += str(i)
                string += "    Mass: "
                string += str(self.getMass(i))
                string += " -> "
                string += str(self.massmap[i])
                print string