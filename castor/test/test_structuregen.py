#! /usr/bin/env python3

# 
# test_structuregen.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import shutil
import tempfile
import os
import math

#MY PATHS

#MY CLASSES
from castor.classes.class_metal import Metal
from castor.classes.class_slab import Slab
from castor.classes.class_vec3 import Vec3_int, Vec3_float
from castor.classes.class_atom import Atom
from castor.classes.class_structure import Structure

#DEFINES
def setup_bulk(metal_in, lc_in, hc_in, packing_in):
    metal = Metal(metal_in, lc_in, hc_in, packing_in)
    dirname = metal.metal + "_" + packing_in + "_Bulk"
    structure = Structure(dirname)
    structure.populate(metal, None)
    return structure

def setup_slab(metal_in, lc_in, hc_in, packing_in, miller_1_in, miller_2_in, miller_3_in, dim_a_in, dim_b_in, layers_in, vacuum_in):
    metal = Metal(metal_in, lc_in, hc_in, packing_in)
    slab = Slab(miller_1_in, miller_2_in, miller_3_in, dim_a_in, dim_b_in, layers_in, vacuum_in)
    dirname = metal.metal + "_" + packing_in + "_" + str(slab.miller.x) + str(slab.miller.y) + str(slab.miller.z) + "_" + str(slab.dim.x) + "x" + str(slab.dim.y) + "_" + str(slab.dim.z) + "L_" + str(slab.vacuum) + "A"
    structure = Structure(dirname)
    structure.populate(metal, slab)
    return structure

def vectortest(vec, check):
    print([vec.x, vec.y, vec.z])
    return all(round(x-y, 12) == 0 for x,y in zip([vec.x, vec.y, vec.z], check))

def test_bulk():
    structure = setup_bulk("Ru", 2.7266, 1.5771, "HCP")
    assert structure.name == "Ru_HCP_Bulk"
    assert structure.lc == 2.7266
    assert vectortest(structure.vec_a,  [1.0, 0.0, 0.0])
    assert vectortest(structure.vec_b,  [0.5, 0.8660254037844387, 0.0])
    assert vectortest(structure.vec_c,  [0.0, 0.0, 1.5771])
    assert structure.elnr == [2]
    assert vectortest(structure.atoms[0].vec, [0.6666666666666666, 0.6666666666666666, 0.75])
    assert vectortest(structure.atoms[1].vec, [0.3333333333333333, 0.3333333333333333, 0.25])

def test_slab():
    structure = setup_slab("CO", 2.4933, 1.6152, "HCP", 0, 0, 1, 3, 3, 5, 15.0)
    assert structure.name == "CO_HCP_001_3x3_5L_15.0A"
    assert structure.lc == 2.4933
    assert vectortest(structure.vec_a,  [1.0, 0.5, 0.0])
    assert vectortest(structure.vec_b,  [-0.5, 0.8660254037844387, 0.0])
    assert vectortest(structure.vec_c,  [0.0, 0.0, 1.6152])
    assert structure.elnr == [2]
    assert vectortest(structure.atoms[0].vec, [0.6666666666666666, 0.6666666666666666, 0.75])
    assert vectortest(structure.atoms[1].vec, [0.3333333333333333, 0.3333333333333333, 0.25])