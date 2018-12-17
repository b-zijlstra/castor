#! /usr/bin/env python

# 
# class_vector.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys, math

#MY PATHS

#MY CLASSES
from .class_rotation import Matrix as rmat

class Vector:
    """Defines a vector of size 3"""
    def __init__(self, xx = 1, yy = 1, zz = 1):
        self.x = float(xx)
        self.y = float(yy)
        self.z = float(zz)
    def __add__(self, rhs):
        return Vector(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    def __sub__(self, rhs):
        return Vector(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    def __mul__(self, rhs):
        if(type(rhs) == int or type(rhs) == float):
            return Vector(self.x * rhs, self.y * rhs, self.z * rhs)
        else:
            return Vector(self.x * rhs.x, self.y * rhs.y, self.z * rhs.z)
    def __div__(self, rhs):
        if(type(rhs) == int or type(rhs) == float):
            temp = 1.0 / rhs
            return Vector(self.x * rhs, self.y * rhs, self.z * rhs)
        else:
            return Vector(self.x / rhs.x, self.y / rhs.y, self.z / rhs.z)
    def length(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
    def length2(self):
        return self.x ** 2 + self.y ** 2 + self.z ** 2
    def triple(self, rhs1, rhs2):
        return (self.x*rhs1.y*rhs2.z +
                self.y*rhs1.z*rhs2.x +
                self.z*rhs1.x*rhs2.y -
                self.z*rhs1.y*rhs2.x -
                self.x*rhs1.z*rhs2.y -
                self.y*rhs1.x*rhs2.z)
    def rotate(self, rmat):
        nx = rmat.g(1,1) * self.x + rmat.g(1,2) * self.y + rmat.g(1,3) * self.z;
        ny = rmat.g(2,1) * self.x + rmat.g(2,2) * self.y + rmat.g(2,3) * self.z;
        nz = rmat.g(3,1) * self.x + rmat.g(3,2) * self.y + rmat.g(3,3) * self.z;
        self.x = nx;
        self.y = ny;
        self.z = nz;
    def getString(self, decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces = spaces_in
        string = ""
        xcor = '{0:.{width}f}'.format(self.x, width=self.decimals)
        ycor = '{0:.{width}f}'.format(self.y, width=self.decimals)
        zcor = '{0:.{width}f}'.format(self.z, width=self.decimals)
        string += '{0:>{width}}'.format(xcor, width=self.decimals+self.spaces+3)
        string += '{0:>{width}}'.format(ycor, width=self.decimals+self.spaces+3)
        string += '{0:>{width}}'.format(zcor, width=self.decimals+self.spaces+3)
        return string
    def write(self, decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces = spaces_in
        string = self.getString(self.decimals, self.spaces)
        print(string)