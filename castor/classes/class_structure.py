#! /usr/bin/env python3

# 
# class_structure.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math

#MY PATHS

#MY CLASSES
from .class_vec3 import Vec3_float
from .class_atom import Atom

class Structure:
    """Defines the poscar content"""
    def __init__(self, name_in):
        self.name = name_in
    def populate(self, metal_in, slab_in):
        if(metal_in.packing=="BCC"):
            self.populate_BCC(metal_in,slab_in)
        elif(metal_in.packing=="FCC"):
            self.populate_FCC(metal_in,slab_in)
        elif(metal_in.packing=="HCP"):
            self.populate_HCP(metal_in,slab_in)
        else:
            print("Unknown packing \"" + metal_in.packing + "\"")
            sys.exit()
    def populate_BCC(self, metal, slab = None):
        self.lc = metal.lc
        if(slab is None):
            self.vec_a = Vec3_float(-0.5, 0.5, 0.5)
            self.vec_b = Vec3_float(0.5, -0.5, 0.5)
            self.vec_c = Vec3_float(0.5, 0.5, -0.5)
            self.elnr = [1]
            self.atoms = [Atom(metal.metal, 0.0, 0.0, 0.0)]
        else:
            print("BCC_Slab method not yet available")
            sys.exit()
    def populate_FCC(self, metal, slab = None):
        self.lc = metal.lc
        if(slab is None):
            self.vec_a = Vec3_float(0.0, 0.5, 0.5)
            self.vec_b = Vec3_float(0.5, 0.0, 0.5)
            self.vec_c = Vec3_float(0.5, 0.5, 0.0)
            self.elnr = [1]
            self.atoms = [Atom(metal.metal, 0.0, 0.0, 0.0)]
        else:
            print("FCC_Slab method not yet available")
            sys.exit()
    def populate_HCP(self, metal, slab = None):
        self.lc = metal.lc
        self.hc = metal.hc
        if(slab is None):
            # self.vec_a = Vec3_float(1.0, 0.5, 0.0)
            # self.vec_b = Vec3_float(-0.5, math.cos(math.pi/6), 0.0)   #use this for 120degree angle at origin
            # self.vec_c = Vec3_float(0.0, 0.0, self.hc)
            self.vec_a = Vec3_float(1.0, 0.0, 0.0)
            self.vec_b = Vec3_float(0.5, math.cos(math.pi/6), 0.0)
            self.vec_c = Vec3_float(0.0, 0.0, self.hc)
            self.elnr = [2]
            # self.atoms = [Atom(metal.metal, 2.0/3.0, 1.0/3.0, 0.25), Atom(metal.metal, 1.0/3.0, 2.0/3.0, 0.75)] #use this atom pos for 120degree angle at origin
            self.atoms = [Atom(metal.metal, 2.0/3.0, 2.0/3.0, 0.75), Atom(metal.metal, 1.0/3.0, 1.0/3.0, 0.25)]
        else:
            self.vec_a = Vec3_float(1.0, 0.5, 0.0)
            self.vec_b = Vec3_float(-0.5, math.cos(math.pi/6), 0.0)
            self.vec_c = Vec3_float(0.0, 0.0, self.hc)
            self.elnr = [2]
            self.atoms = [Atom(metal.metal, 2.0/3.0, 2.0/3.0, 0.75), Atom(metal.metal, 1.0/3.0, 1.0/3.0, 0.25)]


            direct2reciproc_1 = 2.0/(3.0*self.lc*self.lc)*(2.0*slab.miller.x + slab.miller.y)
            direct2reciproc_2 = 2.0/(3.0*self.lc*self.lc)*(slab.miller.x + 2.0*slab.miller.y)
            direct2reciproc_3 = 1.0/(self.lc*self.lc*self.hc*self.hc)*slab.miller.z
            slab.normal.x = self.lc*(direct2reciproc_1*self.vec_a.x + direct2reciproc_2*self.vec_b.x+direct2reciproc_3*self.vec_c.x)
            slab.normal.y = self.lc*(direct2reciproc_1*self.vec_a.y + direct2reciproc_2*self.vec_b.y+direct2reciproc_3*self.vec_c.y)
            slab.normal.z = self.lc*(direct2reciproc_1*self.vec_a.z + direct2reciproc_2*self.vec_b.z+direct2reciproc_3*self.vec_c.z)
            slab.spacing = 1 / math.sqrt(slab.normal.x*slab.normal.x+slab.normal.y*slab.normal.y+slab.normal.z*slab.normal.z)
            
            zerocount = 0
            slab.vec_U = Vec3_float(0.0,0.0,0.0)
            slab.vec_V = Vec3_float(0.0,0.0,0.0)
            if(slab.miller.x==0):
                zerocount += 1
                slab.vec_U = self.vec_a
            if(slab.miller.y==0):
                zerocount += 1
                if(zerocount==1):
                    slab.vec_U = self.vec_b
                else:
                    slab.vec_V = self.vec_b
            if(slab.miller.z==0):
                zerocount += 1
                if(zerocount==1):
                    slab.vec_U = self.vec_c
                else:
                    slab.vec_V = self.vec_c #U and V defined for miller index with 2 zeros

            if(zerocount==3):
                print("Cannot slice 000 surface")
                sys.exit()
            elif(zerocount==1):
                if(slab.miller.x==0):
                    slab.vec_V.x = slab.miller.z*self.vec_b.x - slab.miller.y*self.vec_c.x
                    slab.vec_V.y = slab.miller.z*self.vec_b.y - slab.miller.y*self.vec_c.y
                    slab.vec_V.z = slab.miller.z*self.vec_b.z - slab.miller.y*self.vec_c.z
                if(slab.miller.y==0):
                    slab.vec_V.x = slab.miller.z*self.vec_a.x - slab.miller.x*self.vec_c.x
                    slab.vec_V.y = slab.miller.z*self.vec_a.y - slab.miller.x*self.vec_c.y
                    slab.vec_V.z = slab.miller.z*self.vec_a.z - slab.miller.x*self.vec_c.z
                if(slab.miller.z==0):
                    slab.vec_V.x = slab.miller.y*self.vec_a.x - slab.miller.x*self.vec_b.x
                    slab.vec_V.y = slab.miller.y*self.vec_a.y - slab.miller.x*self.vec_b.y
                    slab.vec_V.z = slab.miller.y*self.vec_a.z - slab.miller.x*self.vec_b.z  #U and V defined for miller index with 1 zero
            elif(zerocount==0):
                slab.vec_U = Vec3_float(1.0*slab.miller.x, 0.5*slab.miller.x, -1.0*self.hc*slab.miller.x)
                slab.vec_V = Vec3_float(-0.5*slab.miller.z, math.cos(math.pi/6)*slab.miller.z, -1.0*self.hc*slab.miller.y) #U and V defined for miller index with no zeros
            
            slab.vec_W = Vec3_float(2.0*slab.dim.z*slab.normal.x, 2.0*slab.dim.z*slab.normal.y, 2.0*slab.dim.z*slab.normal.z)
            # print([slab.vec_U.x, slab.vec_U.y, slab.vec_U.z])
            # print([slab.vec_V.x, slab.vec_V.y, slab.vec_V.z])
            # print([slab.vec_W.x, slab.vec_W.y, slab.vec_W.z])