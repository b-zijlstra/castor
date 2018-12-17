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
def relabel():
    poscar_in = Poscar()
    poscar_out = Poscar()
    poscar_in.read(os.path.dirname(os.path.realpath(__file__)) + "/source/relabel/POSCAR2")
    poscar_out.read(os.path.dirname(os.path.realpath(__file__)) + "/source/relabel/POSCAR")
    poscar_out.relabel(poscar_in, "metal")
    return poscar_out

def vectortest(vec, check):
    print([vec.x, vec.y, vec.z])
    return all(round(x-y, 12) == 0 for x,y in zip([vec.x, vec.y, vec.z], check))

def test_poscar():
    poscar = relabel()
    assert poscar.name == "Co_001_3x3_5L_1x1_15"
    assert vectortest(poscar.unitcell.vec_1,  [0.8660254, 0.5, 0.0])
    assert vectortest(poscar.unitcell.vec_2,  [0.0, 1.0, 0.0])
    assert vectortest(poscar.unitcell.vec_3,  [0.0, 0.0, 3.08218024])
    assert poscar.elnr == [45, 2, 2, 2]
    assert poscar.extra == ['Direct']
    assert vectortest(poscar.atoms[0],  [0.1116411530291996, 0.1094495018142605, 0.4106687751520067])
    assert vectortest(poscar.atoms[21], [0.2249503559256763, 0.5515159503969822, 0.3279375943189227])
    assert vectortest(poscar.atoms[44], [0.8910553131946064, 0.8861626273299941, 0.6729398300143357])
    assert vectortest(poscar.atoms[46], [0.0283259068942601, 0.0802896554572348, 0.2599947017953501])
    assert vectortest(poscar.atoms[48], [0.2095425411536258, 0.0680385721449488, 0.2515655707442391])
    assert vectortest(poscar.atoms[49], [0.9647936000188357, 0.0434528453833594, 0.7785815393652349])
    assert vectortest(poscar.atoms[50], [0.9647936000188357, 0.0434528453833594, 0.2214184606347654])