#! /usr/bin/env python3

# 
# class_frame.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys

#MY PATHS

#MY CLASSES

class Frame:
    """Defines a trajectory frame, which can hold multiple atoms"""
    def __init__(self):
        self.atoms = [] # make an empty list of atoms on initiation
        self.nrat = 0
    def addAtom(self, at):
        self.atoms.append(at)
        self.nrat += 1
    def write(self, comment_in = "Frame", decimals_in = 6, spaces_in = 3):
        print(self.nrat)
        print(comment_in)
        for atom in self.atoms:
            atom.write(decimals_in, spaces_in)