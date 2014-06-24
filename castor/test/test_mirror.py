#! /usr/bin/env python

# 
# test_translate.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 

# #IMPORTS
import sys, os

#MY PATHS

#MY CLASSES
from castor.classes.class_poscar import Poscar

#DEFINES
def mirror_0001():
    poscar_in = Poscar()
    poscar_out = Poscar()
    poscar_in.read(os.path.dirname(os.path.realpath(__file__)) + "/source/mirror/POSCAR")
    poscar_out.read(os.path.dirname(os.path.realpath(__file__)) + "/source/mirror/POSCAR")
    poscar_out.mirror(["ads", "top"], "mirror", "mirror", 0.2)
    poscar_out.move2unitcell()
    return poscar_out

def mirror_1121():
    poscar_in = Poscar()
    poscar_out = Poscar()
    poscar_in.read(os.path.dirname(os.path.realpath(__file__)) + "/source/mirror/CONTCAR")
    poscar_out.read(os.path.dirname(os.path.realpath(__file__)) + "/source/mirror/CONTCAR")
    poscar_out.mirror(["all", "bottom"], "cell_center", "cell_center", 0.2)
    poscar_out.move2unitcell()
    poscar_out.relabel(poscar_in,"metal")
    return poscar_out

def vectortest(vec, check):
    print [vec.x, vec.y, vec.z]
    return all(round(x-y, 12) == 0 for x,y in zip([vec.x, vec.y, vec.z], check))

def test_0001():
    poscar = mirror_0001()
    assert poscar.name == "Co_001_3x3_5L_1x1_15"
    assert vectortest(poscar.unitcell.vec_1,  [0.8660254, 0.5, 0.0])
    assert vectortest(poscar.unitcell.vec_2,  [0.0, 1.0, 0.0])
    assert vectortest(poscar.unitcell.vec_3,  [0.0, 0.0, 3.08218024])
    assert poscar.elnr == [45, 2, 2, 2]
    assert poscar.extra == ['Direct']
    assert vectortest(poscar.atoms[0],  [0.1131629210763434, 0.1093212125207044, 0.4119660507492589])
    assert vectortest(poscar.atoms[21], [0.2138630417052321, 0.5573172891663836, 0.3254910637778055])
    assert vectortest(poscar.atoms[44], [0.8914662503594711, 0.88383223635246, 0.6719051009502915])
    assert vectortest(poscar.atoms[46], [0.3616592402275933, 0.413622988790568, 0.25999470179535056])
    assert vectortest(poscar.atoms[48], [0.5428758744869593, 0.4013719054782823, 0.25156557074423913])
    assert vectortest(poscar.atoms[49], [0.2981269333521692, 0.3767861787166926, 0.7785815393652349])
    assert vectortest(poscar.atoms[50], [0.2981269333521692, 0.3767861787166926, 0.2214184606347651])

def test_1121():
    poscar = mirror_1121()
    assert poscar.name == "Co C O H"
    assert vectortest(poscar.unitcell.vec_1,  [0.81861516, -0.4256520799999999, 0.0])
    assert vectortest(poscar.unitcell.vec_2,  [0.0, 1.0, 0.0])
    assert vectortest(poscar.unitcell.vec_3,  [0.0, 0.0, 2.2293841762382884])
    assert poscar.elnr == [48, 2, 2, 2]
    assert poscar.extra == ['Direct']
    assert vectortest(poscar.atoms[0],  [0.0663476386265634, 0.2558654027454363, 0.338196591576331])
    assert vectortest(poscar.atoms[21], [0.35677928265375436, 0.9711835882095853, 0.6356470412429918])
    assert vectortest(poscar.atoms[44], [0.43067481067023583, 0.24632094586413678, 0.6636999182145391])
    assert vectortest(poscar.atoms[46], [0.9482244221751328, 0.254058349446011, 0.6600538911676725])
    assert vectortest(poscar.atoms[48], [0.7021856095833646, 0.5408239227899432, 0.6644176252969571])
    assert vectortest(poscar.atoms[49], [0.2978143904166354, 0.4591760772100568, 0.3355823747030428])
    assert vectortest(poscar.atoms[50], [0.602304908816041, 0.4619345315882679, 0.7243445633269516])