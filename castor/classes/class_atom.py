#! /usr/bin/env python

# 
# class_atom.py
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
from class_vector import Vector

class Atom:
    """Defines an atom"""
    def __init__(self, ell, xx, yy, zz):
        self.el = ell # the shorthand notation for the element
        self.vec = Vector(xx, yy, zz)
    def translate(self, xx, yy, zz):
        self.vec = self.vec + Vector(xx, yy, zz)
    def place(self, xx, yy, zz):
        self.vec = Vector(xx, yy, zz)
    def scale(self, xx, yy, zz):
        self.vec.x *= xx
        self.vec.y *= yy
        self.vec.z *= zz
    def rotate(self, rmat):
        self.vec.rotate(rmat)
    def setForce(self, fx, fy, fz):
        self.force = Vector(fx, fy, fz)
    def write(self, decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces = spaces_in
        string = '{0:^{width}}'.format(self.el, width=2)
        string += self.vec.getString(self.decimals, self.spaces)
        print string