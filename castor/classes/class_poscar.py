#! /usr/bin/env python

# 
# class_poscar.py
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

#MY PATHS

#MY CLASSES
from class_vector import Vector
# from class_atom import Atom
from class_unitcell import Unitcell

class Poscar:
    """Defines the poscar content"""
    def __init__(self):
        self.name = "NAME"
        self.unitcell = Unitcell()
        self.unitcell.periodic_1 = True
        self.unitcell.periodic_2 = True
        self.unitcell.periodic_3 = True
        self.elnr = []
        self.extra = []
        self.atoms = []
        self.beforeatoms = []
    def read(self, file_in = "POSCAR"):
        with open(file_in, 'r') as inputfile:
            lines = inputfile.readlines()
            linecount = 0
            direct = False
            for line in lines:
                linecount += 1
                if(linecount==1):
                    self.name = line.rstrip()
                    self.beforeatoms.append(line.rstrip())
                elif(linecount==2):
                    match = re.search('^\s*([0-9.-]+)\s*$',line);
                    self.unitcell.lc = float(match.group(1))
                    self.beforeatoms.append(line.rstrip())
                elif(linecount==3):
                    match = re.search('^\s*([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s*$',line);
                    self.unitcell.vec_1 = Vector(float(match.group(1)), float(match.group(2)), float(match.group(3)))
                    self.beforeatoms.append(line.rstrip())
                elif(linecount==4):
                    match = re.search('^\s*([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s*$',line);
                    self.unitcell.vec_2 = Vector(float(match.group(1)), float(match.group(2)), float(match.group(3)))
                    self.beforeatoms.append(line.rstrip())
                elif(linecount==5):
                    match = re.search('^\s*([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s*$',line);
                    self.unitcell.vec_3 = Vector(float(match.group(1)), float(match.group(2)), float(match.group(3)))
                    self.beforeatoms.append(line.rstrip())
                elif(linecount==6):
                    for word in line.split():
                        self.elnr.append(int(word))
                    self.beforeatoms.append(line.rstrip())
                elif(linecount==7):
                    if(re.search('^\s*[Dd].*$',line)):
                        direct = True
                        self.extra.append(line.rstrip())
                    self.beforeatoms.append(line.rstrip())
                elif(linecount==8 and direct==False):
                    if(re.search('^\s*[Dd].*$',line)):
                        direct = True
                        self.extra.append(line.rstrip())
                        self.beforeatoms.append(line.rstrip())
                        continue
                    else:
                        print >> sys.stderr, 'ERROR IN INPUT FILE (not direct?)'
                        sys.exit()
                elif(linecount>7 and direct==True):
                    match = re.search('^\s*([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s*$',line);
                    if(match):
                        self.atoms.append(Vector(float(match.group(1)), float(match.group(2)), float(match.group(3))))
                    else:
                        break
    def write(self, decimals_in = 16, spaces_in = 2, rewrite = False):
        self.decimals = decimals_in
        self.spaces = spaces_in
        if(rewrite==True):
            print self.name
            temp = '{:.{width}f}'.format(self.unitcell.lc, width=self.decimals)
            print '{:>{width}}'.format(temp, width=self.decimals+self.spaces+2)
            print self.unitcell.vec_1.write(self.decimals,self.spaces)
            print self.unitcell.vec_2.write(self.decimals,self.spaces)
            print self.unitcell.vec_3.write(self.decimals,self.spaces)
            firstnum = True
            string = ""
            for number in self.elnr:
                if(firstnum==True):
                    string = str(number)
                    firstnum = False
                else:
                    for i in range(0, self.spaces):
                        string += " "
                    string += str(number)
            print string
            for line in self.extra:
                print line
        else:
            for line in self.beforeatoms:
                print line
        for atom in self.atoms:
            atom.write(self.decimals,self.spaces)
    def translate(self, transX_in, transY_in, sizeX_in, sizeY_in, transZ_in = 0, sizeZ_in = 1):
        for atom in self.atoms:
            atom.x += 1.0*transX_in/sizeX_in
            atom.y += 1.0*transY_in/sizeY_in
            atom.z += 1.0*transZ_in/sizeZ_in
    def rotate(self, rmat):
        for atom in self.atoms:
            vec = self.unitcell.direct2cartesian(Vector(atom.x, atom.y, atom.z))
            vec.rotate(rmat)
            vec = self.unitcell.cartesian2direct(vec)
            atom.x = vec.x
            atom.y = vec.y
            atom.z = vec.z
    def mirror(self, mode_in = ["all", "top"], symmetry_X_in = "mirror", symmetry_Y_in = "mirror", zrange_in = 0.2):
        self.unitcell.invertMatrix()
        z_direct = self.unitcell.cartesian2direct(Vector(0,0,zrange_in)).z
        if(z_direct >= 0.5):
            z_direct = 0.49
        z_max = 0.5 + z_direct
        z_min = 0.5 - z_direct
        # print >> sys.stderr, z_max
        # print >> sys.stderr, z_min

        # cell_center = self.getCellCenter()
        # cell_center.write()
        cell_center = Vector(0.5, 0.5, 0.5)

        atom_center = self.getAtomCenter()
        # atom_center.write()

        sym_x = None
        sym_y = None
        sym_z = 0.5

        if(symmetry_X_in == "mirror"):
            sym_x = None
        elif(symmetry_X_in == "cell_center"):
            sym_x = cell_center.x
        elif(symmetry_X_in == "atom_center"):
            sym_x = atom_center.x
            sym_z = atom_center.z

        if(symmetry_Y_in == "mirror"):
            sym_y = None
        elif(symmetry_Y_in == "cell_center"):
            sym_y = cell_center.y
        elif(symmetry_Y_in == "atom_center"):
            sym_y = atom_center.y
            sym_z = atom_center.z

        newatoms = []
        newelnr = []

        offset = 0
        for el in range(0,len(self.elnr)):
            newel = 0
            if(el == 0 and mode_in[0] == "ads"):
                for i in range(offset, offset+self.elnr[el]):
                    newatoms.append(self.atoms[i])
                    newel += 1
                newelnr.append(newel)
                offset += self.elnr[el]
                continue

            for i in range(offset, offset+self.elnr[el]):
                if((mode_in[1] == "top" and  self.atoms[i].z > z_min) or (mode_in[1] == "bottom" and self.atoms[i].z < z_max)):
                    if(mode_in[1] == "bottom" and self.atoms[i].z < z_min):
                        newatoms.append(self.mirror_atom(self.atoms[i], sym_x, sym_y, sym_z))
                        newel += 1
                    newatoms.append(self.atoms[i])
                    newel += 1
                    if((mode_in[1] == "top" and self.atoms[i].z > z_max)):
                        newatoms.append(self.mirror_atom(self.atoms[i], sym_x, sym_y, sym_z))
                        newel += 1
            newelnr.append(newel)
            offset += self.elnr[el];

        if(newelnr[0] != self.elnr[0]):
            print >> sys.stderr, 'WARNING: NUMBER OF METAL ATOMS CHANGED!'

        self.atoms = newatoms
        self.elnr = newelnr

    def mirror_atom(self, atom_in, sym_x_in = None, sym_y_in = None, sym_z_in = 0.5):
        newatom = Vector(atom_in.x, atom_in.y, atom_in.z)
        if(sym_x_in != None):
            newatom.x = 2.0 * sym_x_in - atom_in.x
        if(sym_y_in != None):
            newatom.y = 2.0 * sym_y_in - atom_in.y
        if(sym_z_in != None):
            newatom.z = 2.0 * sym_z_in - atom_in.z
        return newatom

    def getCellCenter(self):
        vec_a = self.unitcell.direct2cartesian(Vector(1.0, 0.0, 0.0))
        vec_b = self.unitcell.direct2cartesian(Vector(0.0, 1.0, 0.0))
        vec_c = self.unitcell.direct2cartesian(Vector(0.0, 0.0, 1.0))

        center = (vec_a + vec_b + vec_c)*Vector(0.5, 0.5, 0.5)
        cell_center = self.unitcell.cartesian2direct(center)
        return cell_center

    def getAtomCenter(self):
        totvec = Vector(0.0, 0.0, 0.0)
        for atom in self.atoms:
            posvec = self.unitcell.direct2cartesian(Vector(atom.x, atom.y, atom.z))
            totvec += posvec
        div = 1.0/len(self.atoms)
        totvec *= Vector(div, div, div)
        atom_center = self.unitcell.cartesian2direct(totvec)
        return atom_center

    def move2unitcell(self):
        for atom in self.atoms:
            while(atom.x<0):
                atom.x += 1
            while(atom.x>1):
                atom.x -= 1
            while(atom.y<0):
                atom.y += 1
            while(atom.y>1):
                atom.y -= 1
            while(atom.z<0):
                atom.z += 1
            while(atom.z>1):
                atom.z -= 1
    def relabel(self,poscar_in, target = "metal"):
        if((target == "all" and self.elnr!=poscar_in.elnr) or (target == "metal" and self.elnr[0]!=poscar_in.elnr[0]) or (target == "adsorbate" and self.elnr[1:]!=poscar_in.elnr[1:])):
            print >> sys.stderr, 'Number of elements do not match for relabeling'
            sys.exit()

        dictionary1 = {}
        dictionary2 = {}
        offset = 0
        difference_max = self.unitcell.direct2cartesian(Vector(0.5,0.5,0.5)).length()
        for el in range(0,len(self.elnr)):
            for i in range(offset, offset+self.elnr[el]):
                difference = difference_max
                foundmatch = False
                if(target=="all" or (target == "metal" and el==0) or (target == "adsorbate" and el>0)):
                    for j in range(offset, offset+self.elnr[el]):
                        if((j in dictionary2) == False):
                            diff = self.atoms[i]-poscar_in.atoms[j]
                            diff = self.unitcell.getModularDirect(diff)
                            diff = self.unitcell.direct2cartesian(diff)
                            diff = diff.length()
                            if(diff<difference):
                                difference = diff
                                match = j
                                foundmatch = True
                if(foundmatch==True):
                    # print >> sys.stderr, str(i) +  "--" + str(match)
                    dictionary1[i]= match
                    dictionary2[match]= i
                else:
                    if((i in dictionary2) == False):
                        # print >> sys.stderr, str(i) +  "--" + str(i)
                        dictionary1[i]= i
                        dictionary2[i]= i
                    else:
                        print >> sys.stderr, 'Failed to create label dictionary'
                        sys.exit()

            offset += self.elnr[el];

        tempatoms = []
        for i in range(0, len(self.atoms)):
            tempatoms.append(self.atoms[dictionary2[i]])
        self.atoms = tempatoms