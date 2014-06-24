#! /usr/bin/env python

# 
# class_unitcell.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import numpy as np

#MY PATHS

#MY CLASSES
from class_vector import Vector


class Unitcell:
    """Defines a Unitcell"""
    def __init__(self, vec1 = Vector(100, 0, 0), vec2 = Vector(0, 100, 0), vec3 = Vector(0, 0, 100), lc_in = 1.0):
        self.lc = lc_in
        self.vec_1 = vec1
        self.vec_2 = vec2
        self.vec_3 = vec3
        self.periodic_1 = False
        self.periodic_2 = False
        self.periodic_3 = False
    def invertMatrix(self):
        self.vec_1_lc = self.vec_1 * self.lc
        self.vec_2_lc = self.vec_2 * self.lc
        self.vec_3_lc = self.vec_3 * self.lc
        unitcelmatrix = [ [self.vec_1_lc.x, self.vec_1_lc.y, self.vec_1_lc.z], [self.vec_2_lc.x, self.vec_2_lc.y, self.vec_2_lc.z], [self.vec_3_lc.x, self.vec_3_lc.y, self.vec_3_lc.z] ]
        unitcelmatrix_inv = np.linalg.inv(unitcelmatrix)
        self.vec_1_inv = Vector(unitcelmatrix_inv[0][0], unitcelmatrix_inv[0][1], unitcelmatrix_inv[0][2])
        self.vec_2_inv = Vector(unitcelmatrix_inv[1][0], unitcelmatrix_inv[1][1], unitcelmatrix_inv[1][2])
        self.vec_3_inv = Vector(unitcelmatrix_inv[2][0], unitcelmatrix_inv[2][1], unitcelmatrix_inv[2][2])
    def direct2cartesian(self, vec_in):
        direct = vec_in
        cartesian = Vector()
        cartesian.x = self.lc * (direct.x * self.vec_1.x + direct.y * self.vec_2.x + direct.z * self.vec_3.x)
        cartesian.y = self.lc * (direct.x * self.vec_1.y + direct.y * self.vec_2.y + direct.z * self.vec_3.y)
        cartesian.z = self.lc * (direct.x * self.vec_1.z + direct.y * self.vec_2.z + direct.z * self.vec_3.z)
        return cartesian
    def cartesian2direct(self, vec_in, reinvert=False):
        cartesian = vec_in
        direct = Vector()
        if(reinvert==True):
            self.invertMatrix()
        direct.x = cartesian.x * self.vec_1_inv.x + cartesian.y * self.vec_2_inv.x + cartesian.z * self.vec_3_inv.x
        direct.y = cartesian.x * self.vec_1_inv.y + cartesian.y * self.vec_2_inv.y + cartesian.z * self.vec_3_inv.y
        direct.z = cartesian.x * self.vec_1_inv.z + cartesian.y * self.vec_2_inv.z + cartesian.z * self.vec_3_inv.z
        return direct;
    def getModularDirect(self, vec_in, mod_in = Vector(1.0, 1.0, 1.0)):
        if(self.periodic_1 == True or self.periodic_2 == True or self.periodic_3 == True):
            direct = vec_in
            modvec = mod_in
            result = Vector(direct.x, direct.y, direct.z)
            if(self.periodic_1 == True):
                result.x -= modvec.x * round(direct.x/modvec.x)
            if(self.periodic_2 == True):
                result.y -= modvec.y * round(direct.y/modvec.y)
            if(self.periodic_3 == True):
                result.z -= modvec.z * round(direct.z/modvec.z)
            return result
        else:
            return vec_in
    def getModularCartesian(self, vec_in):
        if(self.periodic_1 == True or self.periodic_2 == True or self.periodic_3 == True):
            cartesian = vec_in
            direct = self.cartesian2direct(cartesian)
            direct = self.getModularDirect(direct)
            cartesian = self.direct2cartesian(direct)
            return cartesian
        else:
            return vec_in
    def getVolume(self):
        return abs(self.vec_1.triple(self.vec_2,self.vec_3)) * self.lc * self.lc * self.lc
    def write(self, decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces = spaces_in
        self.vec_1_lc = self.vec_3 * self.lc
        self.vec_2_lc = self.vec_3 * self.lc
        self.vec_3_lc = self.vec_3 * self.lc
        string = self.vec_1_lc.getString(self.decimals, self.spaces)
        if(self.periodic_1 == True):
            string += "   (periodic)"
        print string
        string = self.vec_2_lc.getString(self.decimals, self.spaces)
        if(self.periodic_2 == True):
            string += "   (periodic)"
        print string
        string = self.vec_3_lc.getString(self.decimals, self.spaces)
        if(self.periodic_3 == True):
            string += "   (periodic)"
        print string