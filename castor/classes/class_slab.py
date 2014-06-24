#! /usr/bin/env python

# 
# class_slab.py
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
from class_vec3 import Vec3_int, Vec3_float

class Slab:
    """Defines the surface slab"""
    def __init__(self, miller_1_in, miller_2_in, miller_3_in, dim_a_in, dim_b_in, layers_in, vacuum_in):
        self.miller = Vec3_int(int(miller_1_in), int(miller_2_in), int(miller_3_in))
        self.dim = Vec3_int(int(dim_a_in), int(dim_b_in), int(layers_in))
        self.vacuum = float(vacuum_in)
        self.normal = Vec3_float(0.0, 0.0, 0.0)
        self.spacing = 0.0