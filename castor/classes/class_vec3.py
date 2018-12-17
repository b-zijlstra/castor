#! /usr/bin/env python

# 
# class_vec3.py
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

class Vec3_int:
    """Defines vector of 3 integers"""
    def __init__(self, int1_in, int2_in, int3_in):
        self.x = int(int1_in)
        self.y = int(int2_in)
        self.z = int(int3_in)
    def scale(self, xx, yy = None, zz = None):
        if((yy == None and zz != None) or (yy != None and zz == None)):
            print('scale funtion of Vec3_int has wrong amount of arguments', file=sys.stderr)
            sys.exit()
        elif(yy is None and zz is None):
            self.x *= int(xx)
            self.y *= int(xx)
            self.z *= int(xx)
        else:
            self.x *= int(xx)
            self.y *= int(yy)
            self.z *= int(zz)
    def add(self, xx, yy = None, zz = None):
        if((yy == None and zz != None) or (yy != None and zz == None)):
            print('add funtion of Vec3_int has wrong amount of arguments', file=sys.stderr)
            sys.exit()
        elif(yy is None and zz is None):
            self.x += int(xx)
            self.y += int(xx)
            self.z += int(xx)
        else:
            self.x += int(xx)
            self.y += int(yy)
            self.z += int(zz)
    def addVec(self, vec_in):
        self.x += int(vec_in.x)
        self.y += int(vec_in.y)
        self.z += int(vec_in.z)

class Vec3_float:
    """Defines vector of 3 floats"""
    def __init__(self, float1_in, float2_in, float3_in):
        self.x = float(float1_in)
        self.y = float(float2_in)
        self.z = float(float3_in)
    def write(self, decimals_in, spaces_in):
        string = ""
        for number in [self.x, self.y, self.z]:
            temp = '{0:.{width}f}'.format(number, width=decimals_in)
            string += '{0:>{width}}'.format(temp, width=decimals_in+spaces_in+2)
        return string
    def scale(self, xx, yy = None, zz = None):
        if((yy == None and zz != None) or (yy != None and zz == None)):
            print('scale funtion of Vec3_float has wrong amount of arguments', file=sys.stderr)
            sys.exit()
        elif(yy is None and zz is None):
            self.x *= float(xx)
            self.y *= float(xx)
            self.z *= float(xx)
        else:
            self.x *= float(xx)
            self.y *= float(yy)
            self.z *= float(zz)
    def add(self, xx, yy = None, zz = None):
        if((yy == None and zz != None) or (yy != None and zz == None)):
            print('add funtion of Vec3_float has wrong amount of arguments', file=sys.stderr)
            sys.exit()
        elif(yy is None and zz is None):
            self.x += float(xx)
            self.y += float(xx)
            self.z += float(xx)
        else:
            self.x += float(xx)
            self.y += float(yy)
            self.z += float(zz)
    def addVec(self, vec_in):
        self.x += float(vec_in.x)
        self.y += float(vec_in.y)
        self.z += float(vec_in.z)