#! /usr/bin/env python

# 
# class_trajectory.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import re

#MY PATHS

#MY CLASSES
from .class_unitcell import Unitcell
from .class_frame import Frame
from .class_atom import Atom
from .class_vector import Vector

class Trajectory:
    """Defines a trajectory, which can hold multiple frames and a unitcell"""
    def __init__(self):
        self.frames = [] # make an empty list of frames on initiation
        self.xyz = None
        self.unitcell = Unitcell()
        self.coordinates = "cartesian"

    def readXYZ(self, xyz_in):
        with open(xyz_in, 'r') as inputfile:
            self.frame = Frame()
            matrixlist_frames = []
            stringlist_elements = []
            self.atomcount = 0
            self.numberofatoms = 0
            self.readmode = "empty"
            linecount = 0
            for line in inputfile:
                linecount += 1
                if(self.readmode=="empty"):
                    if(self.checkEmpty(line) == True):
                        continue
                    number = self.checkNumberAtoms(line)
                    if(number != None):
                        self.atomcount = 0
                        self.frame = Frame()
                        self.readmode = "number"
                        self.numberofatoms = number
                        continue
                    if(self.checkVectorXYZ(line) == True):
                        continue
                    atom = self.checkAtomXYZ(line)
                    if(atom == None):
                        print("WARNING - No method available for line " + str(linecount) + ":")
                        print("--- " + line)
                        continue
                    else:
                        print("WARNING - Atom entry on line " + str(linecount) + " does not belong to any structure")
                        print("--- " + line)
                        continue

                if(self.readmode=="number"):
                    if(self.checkEmpty(line) == True):
                        continue
                    atom = self.checkAtomXYZ(line)
                    if(atom != None):
                        self.atomcount += 1
                        self.readmode = "xyz"
                        self.frame.addAtom(atom)
                        self.checkEndOfFrame()
                        continue

                if(self.readmode=="xyz"):
                    if(self.checkEmpty(line) == True):
                        continue
                    atom = self.checkAtomXYZ(line)
                    if(atom != None):
                        self.atomcount += 1
                        self.frame.addAtom(atom)
                        self.checkEndOfFrame()
                        continue
                    number = checkNumberAtoms(line)
                    if(number != None):
                        print("WARNING - Adsorbate structure " + str(len(self.self.frames)+1) + " has less atoms than defined. Ending prematurely!")
                        self.frames.append(self.frame)
                        self.frame = Frame()
                        self.atomcount = 0
                        self.readmode = "number"
                        self.numberofatoms = number
                        continue
                    else:
                        print("WARNING - No method available for line " + str(linecount) + ":")
                        print("--- " + line)
                        continue

            if(self.readmode != "empty"):
                print("WARNING - Adsorbate structure " + str(len(self.frames)+1) + " has less atoms than defined. Ending prematurely at end of file!")
                self.frames.append(self.frame)

    def checkEmpty(self, string_in):
        match = re.search('^[ \t]*$', string_in)
        if(match):
            return True
        else:
            return False;

    def checkNumberAtoms(self, string_in):
        match = re.search('^[ \t]*([0-9]+)[ \t]*$',string_in)
        if(match):
            return int(match.group(1))
        else:
            return None;

    def checkAtomXYZ(self, string_in):
        match = re.search('^[ \t]*([A-Za-z]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+).*$',string_in)
        if(match):
            el = match.group(1)
            x = match.group(2)
            y = match.group(3)
            z = match.group(4)
            atom = Atom(el,x,y,z)
            return atom
        else:
            return None

    def checkVectorXYZ(self, string_in):
        match = re.search('^[ \t]*VEC1[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+).*$',string_in)
        if(match):
            self.unitcell.vec_1 = Vector(match.group(1),match.group(2),match.group(3))
            return True
        match = re.search('^[ \t]*VEC2[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+).*$',string_in)
        if(match):
            self.unitcell.vec_2 = Vector(match.group(1),match.group(2),match.group(3))
            return True
        match = re.search('^[ \t]*VEC3[ \t]+([0-9.-]+)[ \t]+([0-9.-]+)[ \t]+([0-9.-]+).*$',string_in)
        if(match):
            self.unitcell.vec_3 = Vector(match.group(1),match.group(2),match.group(3))
            return True
        return False

    def checkEndOfFrame(self):
        if(self.atomcount >= self.numberofatoms):
            self.frames.append(self.frame)
            self.frame = Frame()
            self.atomcount = 0
            self.readmode = "empty"
            self.numberofatoms = 0

    def setBoundaries(self, bound1, bound2, bound3):
        self.unitcell.periodic_1 = bool(bound1)
        self.unitcell.periodic_2 = bool(bound2)
        self.unitcell.periodic_3 = bool(bound3)
        if(self.unitcell.periodic_1 == True or self.unitcell.periodic_2 == True or self.unitcell.periodic_3 == True):
            self.unitcell.invertMatrix()

    def write(self, printmode = "all", decimals_in = 6, spaces_in = 3):
        if(printmode == "all"):
            count = 0
            for frame in self.frames:
                count += 1
                comment = "Frame "+str(count)
                frame.write(comment, decimals_in, spaces_in)
        else:
            print(str(len(self.frames)) + " frames read from " + self.xyz)

    def getLastElnames(self):
        frame = self.frames[-1]
        atomnamelist = []
        atomname = None
        for atom in frame.atoms:
            if(atom.el == atomname):
                pass
            else:
                atomnamelist.append(atom.el)
                atomname = atom.el
        return atomnamelist

    def getLastElnr(self):
        frame = self.frames[-1]
        atomnrlist = []
        atomname = None
        atomcount = 0
        for atom in frame.atoms:
            if(atomname == None):
                atomname = atom.el
                atomcount = 1
            elif(atom.el == atomname):
                atomcount += 1
            else:
                atomnrlist.append(atomcount)
                atomname = atom.el
                atomcount = 1
        atomnrlist.append(atomcount)
        return atomnrlist

    def getLastCoord(self):
        frame = self.frames[-1]
        atomvectors = []
        for atom in frame.atoms:
            atomvectors.append(atom.vec)
        return atomvectors