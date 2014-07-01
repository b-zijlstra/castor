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

class Hessian:
    """Defines a Hessian"""
    def __init__(self, temp_in = 800):
        self.decimals = 6
        self.spaces = 3
        self.vector_atoms = np.dtype(int)
        self.vector_axes = np.dtype(str)
        self.vector_names = np.dtype(str)
        self.vector_elnr = np.dtype(int)
        self.vector_elements = np.dtype(str)
        self.vector_masses = np.dtype(float)
        self.matrix_nonsym = None
        self.matrix_sym = None
        self.matrix_mass_ori = None
        self.matrix_mass_new = None
        self.matrix_diag_ori = None
        self.matrix_diag_new = None
        self.vector_THz_ori = np.dtype(float)
        self.vector_THz_new = np.dtype(float)
        self.vector_2PiTHz_ori = np.dtype(float)
        self.vector_2PiTHz_new = np.dtype(float)
        self.vector_cm1_ori = np.dtype(float)
        self.vector_cm1_new = np.dtype(float)
        self.vector_meV_ori = np.dtype(float)
        self.vector_meV_new = np.dtype(float)
        self.vector_imaginary_ori = np.dtype(bool)
        self.vector_imaginary_new = np.dtype(bool)
        self.numberset = set()
        self.massmap = dict()
        self.outcar = None
        self.numbers = None
        self.element = None
        self.mass = None
        self.changes = False
        self.temp = temp_in # in Kelvin - Used for the vibrational partition function
        self.kb = 8.617332478e-2 # Botzmann constant in meV/K
    def read(self, outcar_in):
        self.outcar = outcar_in
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
            matrix_tot = np.matrix(floatlist_list_rows)
            matrix_size = len(matrix_tot)/3
            self.matrix_nonsym = matrix_tot[0 : matrix_size]
            self.matrix_sym = matrix_tot[matrix_size : matrix_size*2]
            self.matrix_mass_ori = matrix_tot[matrix_size*2 : matrix_size*3]
            self.vector_atoms = intlist_atoms[0:matrix_size]
            self.vector_axes = charlist_axes[0:matrix_size]
            self.vector_elements = stringlist_elements
            self.vector_masses = floatlist_masses
        if(self.matrix_sym == None or len(self.matrix_sym) == 0 or self.matrix_sym[0].size == 0):
            print "Reading from " + outcar_in + " failed! Could not find hessian matrix."
            sys.exit()
        if(self.matrix_mass_ori == None or len(self.matrix_mass_ori) == 0 or self.matrix_mass_ori[0].size == 0):
            print "Reading from " + outcar_in + " failed! Could not find hessian matrix."
            sys.exit()
        if(self.vector_elnr == None):
            print "Reading from " + outcar_in + " failed! Could not find number of elements."
            sys.exit()

    def mapMass(self, element_in = None, mass_in = None, numbers_in = None):
        self.numbers = numbers_in
        self.element = element_in
        self.mass = mass_in

        self.numberset = self.string2numberset(numbers_in, element_in)
        
        self.massmap.clear()
        for i in range(0,len(self.vector_atoms)):
            if(self.vector_atoms[i] not in self.massmap):
                if(self.vector_atoms[i] in self.numberset and mass_in != None):
                    if(self.changes == False and mass_in != self.getMass(self.vector_atoms[i])):
                        self.changes = True
                    self.massmap[self.vector_atoms[i]] = mass_in
                else:
                    self.massmap[self.vector_atoms[i]] = self.getMass(self.vector_atoms[i])

    def newMassMatrix(self):
        matrix_mass_new = self.matrix_sym

        for i in range(0,len(matrix_mass_new)):
            mass_i = self.massmap[self.vector_atoms[i]]
            for j in range(0,matrix_mass_new[i].size):
                mass_j = self.massmap[self.vector_atoms[j]]
                matrix_mass_new[i,j] /= math.sqrt(mass_i*mass_j)
        self.matrix_mass_new = matrix_mass_new 

    def string2numberset(self, string_in, element_in):
        numberset = set()
        if(string_in == None or string_in == "all"):
            for x in range(0,len(self.vector_atoms)):
                if(self.getElement(self.vector_atoms[x])==element_in):
                    numberset.add(self.vector_atoms[x])
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
                                x = range(int(x[0]),np.max(self.vector_atoms)+1)
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
        for i in range(0,len(self.vector_elnr)):
            atomsum += self.vector_elnr[i]
            if(number_in <= atomsum):
                return self.vector_elements[i]
        print "Could not get element name"
        sys.exit()

    def getMass(self, number_in):
        atomsum = 0
        for i in range(0,len(self.vector_elnr)):
            atomsum += self.vector_elnr[i]
            if(number_in <= atomsum):
                return self.vector_masses[i]
        print "Could not get element mass"
        sys.exit()

    def diagonalize(self):
        self.diagonalize_ori()
        self.diagonalize_new()

    def diagonalize_ori(self):
        if(self.matrix_mass_ori == None or len(self.matrix_mass_ori) == 0 or self.matrix_mass_ori[0].size == 0):
            print "Can not diagonalize original hessian matrix because it does not exist!"
            sys.exit()
        else:
            self.matrix_diag_ori = np.linalg.eigh(self.matrix_mass_ori)
    
    def diagonalize_new(self):
        if(self.matrix_mass_new == None or len(self.matrix_mass_new) == 0 or self.matrix_mass_new[0].size == 0):
            print "Can not diagonalize original hessian matrix because it does not exist!"
            sys.exit()
        else:
            self.matrix_diag_new = np.linalg.eigh(self.matrix_mass_new)

    def calcFreqs(self):
        self.calcFreqs_ori()
        self.calcFreqs_new()

    def calcFreqs_ori(self):
        eig_to_THz = 15.633304592
        eig_to_2PiTHz = 98.2269497148
        eig_to_cm1 = 521.47091
        eig_to_meV = 64.6541499
        list_THz_ori = []
        list_2PiTHz_ori = []
        list_cm1_ori = []
        list_meV_ori = []
        list_imaginary_ori = []
        for i in range(0,len(self.matrix_diag_ori[0])):
            eigenval = self.matrix_diag_ori[0][i]
            if(eigenval < 0):
                list_imaginary_ori.append(False)
                eigenval = math.sqrt(abs(eigenval))
            else:
                list_imaginary_ori.append(True)
                eigenval = math.sqrt(eigenval)
            list_THz_ori.append(eigenval*eig_to_THz)
            list_2PiTHz_ori.append(eigenval*eig_to_2PiTHz)
            list_cm1_ori.append(eigenval*eig_to_cm1)
            list_meV_ori.append(eigenval*eig_to_meV)
        self.vector_THz_ori = list_THz_ori
        self.vector_2PiTHz_ori = list_2PiTHz_ori
        self.vector_cm1_ori = list_cm1_ori
        self.vector_meV_ori = list_meV_ori
        self.vector_imaginary_ori = list_imaginary_ori

    def calcFreqs_new(self):
        eig_to_THz = 15.633304592
        eig_to_2PiTHz = 98.2269497148
        eig_to_cm1 = 521.47091
        eig_to_meV = 64.6541499
        list_THz_new = []
        list_2PiTHz_new = []
        list_cm1_new = []
        list_meV_new = []
        list_imaginary_new = []
        for i in range(0,len(self.matrix_diag_new[0])):
            eigenval = self.matrix_diag_new[0][i]
            if(eigenval < 0):
                list_imaginary_new.append(False)
                eigenval = math.sqrt(abs(eigenval))
            else:
                list_imaginary_new.append(True)
                eigenval = math.sqrt(eigenval)
            list_THz_new.append(eigenval*eig_to_THz)
            list_2PiTHz_new.append(eigenval*eig_to_2PiTHz)
            list_cm1_new.append(eigenval*eig_to_cm1)
            list_meV_new.append(eigenval*eig_to_meV)
        self.vector_THz_new = list_THz_new
        self.vector_2PiTHz_new = list_2PiTHz_new
        self.vector_cm1_new = list_cm1_new
        self.vector_meV_new = list_meV_new
        self.vector_imaginary_new = list_imaginary_new

    def write(self, printmode):
        kbT = self.temp * self.kb
        if(printmode == "all"):
            print "---------------------------"
            print "-        Settings:        -"
            print "---------------------------"
            print "Outcar used = " + str(self.outcar)
            print "Element to set new mass = " + str(self.element)
            if(self.numbers == None):
                print "Atom numbers of element to set new mass = all"
            else:
                print "Atom numbers of element to set new mass = " + str(self.numbers)
            print "New mass = " + str(self.mass)
            print "\n"
            if(self.matrix_mass_ori != None):
                print "---------------------------"
                print "-  Original mass matrix:  -"
                print "---------------------------"
                self.printMatrix(self.matrix_mass_ori)
                print "\n"
            if(self.matrix_mass_new != None):
                print "---------------------------"
                print "-        Changes:         -"
                print "---------------------------"
                self.printChanges()
                print "\n"
                print "---------------------------"
                print "-    New mass matrix:     -"
                print "---------------------------"
                self.printMatrix(self.matrix_mass_new)
                print "\n"
            if(len(self.vector_meV_ori) != 0):
                print "---------------------------"
                print "-  Original frequencies:  -"
                print "---------------------------"
                for i in range(0,len(self.vector_meV_ori)):
                    string = ""
                    string += '{:>{width}}'.format(i+1, width=self.spaces+1)
                    if(self.vector_imaginary_ori[i]==False):
                        string += " f  ="
                    elif(self.vector_imaginary_ori[i]==True):
                        string += " f/i="
                    data = '{:.{width}f} THz'.format(self.vector_THz_ori[i], width=self.decimals)
                    string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+6)
                    data = '{:.{width}f} 2PiTHz'.format(self.vector_2PiTHz_ori[i], width=self.decimals)
                    string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+12)
                    data = '{:.{width}f} cm-1'.format(self.vector_cm1_ori[i], width=self.decimals)
                    string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+11)
                    data = '{:.{width}f} meV'.format(self.vector_meV_ori[i], width=self.decimals)
                    string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+9)
                    print string
                print "\n"
                zpe = 0.0
                nu = 1.0
                for i in range(0,len(self.vector_meV_ori)):
                    if(self.vector_imaginary_ori[i]==False):
                        zpe += self.vector_meV_ori[i]
                        nu *= 1.0 / (1.0 - math.exp(-self.vector_meV_ori[i] / kbT))
                print 'Total ZPE contribution of frequencies in eV: {:.{width}f}'.format(zpe/2000, width=self.decimals)
                print 'Vibrational partition function: {:.{width}f}'.format(nu, width=self.decimals)
                print "\n"
            if(len(self.vector_meV_new) != 0):
                print "---------------------------"
                print "-    New frequencies:     -"
                print "---------------------------"
                for i in range(0,len(self.vector_meV_new)):
                    string = ""
                    string += '{:>{width}}'.format(i+1, width=self.spaces+1)
                    if(self.vector_imaginary_new[i]==False):
                        string += " f  ="
                    elif(self.vector_imaginary_new[i]==True):
                        string += " f/i="
                    data = '{:.{width}f} THz'.format(self.vector_THz_new[i], width=self.decimals)
                    string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+6)
                    data = '{:.{width}f} 2PiTHz'.format(self.vector_2PiTHz_new[i], width=self.decimals)
                    string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+12)
                    data = '{:.{width}f} cm-1'.format(self.vector_cm1_new[i], width=self.decimals)
                    string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+11)
                    data = '{:.{width}f} meV'.format(self.vector_meV_new[i], width=self.decimals)
                    string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+9)
                    print string
                print "\n"
                zpe = 0.0
                nu = 1.0
                for i in range(0,len(self.vector_meV_new)):
                    if(self.vector_imaginary_new[i]==False):
                        zpe += self.vector_meV_new[i]
                        nu *= 1.0 / (1.0 - math.exp(-self.vector_meV_new[i] / kbT))
                print 'Total ZPE contribution of frequencies in eV: {:.{width}f}'.format(zpe/2000, width=self.decimals)
                print 'Vibrational partition function: {:.{width}f}'.format(nu, width=self.decimals)
                print "\n"
        elif(len(self.vector_meV_ori) != 0):
            # print "-----------"
            zpe = 0.0
            nu = 1.0
            for i in range(0,len(self.vector_meV_ori)):
                if(self.vector_imaginary_ori[i]==False):
                    zpe += self.vector_meV_ori[i]
                    nu *= 1.0 / (1.0 - math.exp(-self.vector_meV_ori[i] / kbT))
            print 'Total ZPE contribution of frequencies in eV: {:.{width}f}'.format(zpe/2000, width=self.decimals)
            print 'Vibrational partition function: {:.{width}f}'.format(nu, width=self.decimals)

            if(len(self.vector_meV_new) != 0):
                # print "---"
                self.printChanges()
                # print "---"
                zpe = 0.0
                nu = 1.0
                for i in range(0,len(self.vector_meV_new)):
                    if(self.vector_imaginary_new[i]==False):
                        zpe += self.vector_meV_new[i]
                        nu *= 1.0 / (1.0 - math.exp(-self.vector_meV_new[i] / kbT))
                print 'Total ZPE contribution of new frequencies in eV: {:.{width}f}'.format(zpe/2000, width=self.decimals)
                print 'New vibrational partition function: {:.{width}f}'.format(nu, width=self.decimals)
            # print "-----------"

    def printMatrix(self, matrix_in):
        self.printAxis_X()
        for i in range(0,len(matrix_in)):
            name = str(self.vector_atoms[i])+self.vector_axes[i]
            string = '{:^{width}}'.format(name, width=self.spaces+2)
            for j in range(0,matrix_in[i].size):
                data = '{:.{width}f}'.format(matrix_in[i,j], width=self.decimals)
                string += '{:>{width}}'.format(data, width=self.decimals+self.spaces+3)
            print string

    def printAxis_X(self):
        string = '{:^{width}}'.format("", width=self.spaces+2)
        for i in range(0,len(self.vector_atoms)):
            name = str(self.vector_atoms[i])+self.vector_axes[i]
            string += '{:>{width}}'.format(name, width=self.decimals+self.spaces+3)
        print string

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