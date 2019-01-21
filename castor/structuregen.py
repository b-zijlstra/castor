#! /usr/bin/env python3

# 
# structuregen.py
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
from classes.class_metal import Metal
from classes.class_slab import Slab
from classes.class_vec3 import Vec3_int, Vec3_float
from classes.class_atom import Atom
from classes.class_folder import Folder
from classes.class_formatting import Formatting
from classes.class_structure import Structure

#DEFINES
def setup_bulk(metal_in, lc_in, hc_in, packing_in):
    metal = Metal(metal_in, lc_in, hc_in, packing_in)
    dirname = metal.metal + "_" + packing_in + "_Bulk"
    formatting = Formatting(16,4)
    structure = Structure(dirname)
    structure.populate(metal, None)
    make_folder(dirname, structure, formatting)
    print("DONE!")

def setup_slab(metal_in, lc_in, hc_in, packing_in, miller_1_in, miller_2_in, miller_3_in, dim_a_in, dim_b_in, layers_in, vacuum_in):
    metal = Metal(metal_in, lc_in, hc_in, packing_in)
    slab = Slab(miller_1_in, miller_2_in, miller_3_in, dim_a_in, dim_b_in, layers_in, vacuum_in)
    dirname = metal.metal + "_" + packing_in + "_" + str(slab.miller.x) + str(slab.miller.y) + str(slab.miller.z) + "_" + str(slab.dim.x) + "x" + str(slab.dim.y) + "_" + str(slab.dim.z) + "L_" + str(slab.vacuum) + "A"
    formatting = Formatting(16,4)
    structure = Structure(dirname)
    structure.populate(metal, slab)
    make_folder(dirname, structure, formatting)
    print("DONE!")

def populate(metal_in, slab_in):
    if(metal_in.packing=="BCC"):
        structure.populate_BCC(metal_in,slab_in)
    elif(metal_in.packing=="FCC"):
        structure.populate_FCC(metal_in,slab_in)
    elif(metal_in.packing=="HCP"):
        structure.populate_HCP(metal_in,slab_in)
    else:
        print("Unknown packing \"" + metal_in.packing + "\"")
        sys.exit()

def make_folder(dirname_in, structure_in, formatting_in):
    folder = Folder(dirname_in)
    folder.make_root()
    folder.add_file("source/INCAR", "INCAR")
    folder.add_file("source/KPOINTS" , "KPOINTS")
    folder.add_file("source/POTCAR", "POTCAR")
    folder.add_qsub("source/submit.qsub", "submit.qsub", dirname_in)
    folder.add_poscar(structure_in, formatting_in, "POSCAR")



#EXECUTION
if len(sys.argv)==5:
    arguments = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    setup_bulk(*arguments)
elif len(sys.argv)==12:
    arguments = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]), int(sys.argv[10]), float(sys.argv[11])
    setup_slab(*arguments)
else:
    print("Unexpected amount of arguments.")