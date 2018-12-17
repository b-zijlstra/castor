#! /usr/bin/env python

# 
# test_rotate.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys, os

#MY PATHS

#MY CLASSES
from castor.classes.class_poscar import Poscar
from castor.classes.class_rotation import Matrix, AxisAngle
from castor.classes.class_vector import Vector

#DEFINES
def rotate(rotX1_in, rotY1_in, rotZ1_in, rotX2_in, rotY2_in, rotZ2_in, angle_in, input_in):
    metalOnly = True
    poscar = Poscar()
    poscar.read(input_in)
    poscar.unitcell.invertMatrix()
    vx = rotX2_in - rotX1_in
    vy = rotY2_in - rotY1_in
    vz = rotZ2_in - rotZ1_in
    newvec = poscar.unitcell.direct2cartesian(Vector(vx,vy,vz))
    q = AxisAngle(newvec.x, newvec.y, newvec.z, angle_in)
    poscar.translate(-rotX2_in, -rotY2_in, 1, 1, -rotZ2_in, 1)
    poscar.rotate(q.toMatrix())
    poscar.translate(rotX2_in, rotY2_in, 1, 1, rotZ2_in, 1)
    return poscar

def vectortest(vec, check):
    print([vec.x, vec.y, vec.z])
    return all(round(x-y, 12) == 0 for x,y in zip([vec.x, vec.y, vec.z], check))

def test_poscar():
    poscar = rotate(0.0, 0.5, 0.5, 1.0, 0.5, 0.5, 45, os.path.dirname(os.path.realpath(__file__)) + "/source/rotate/POSCAR")
    assert poscar.name == "polythiophene"
    assert vectortest(poscar.unitcell.vec_1, [9.0000000000000000, 0.0000000000000000, 0.0000000000000000])
    assert vectortest(poscar.unitcell.vec_2, [0.0000000000000000, 8.0000000000000000, 0.0000000000000000])
    assert vectortest(poscar.unitcell.vec_3, [0.0000000000000000, 0.0000000000000000, 8.0000000000000000])
    assert poscar.elnr == [8, 4, 2]
    assert poscar.extra == ['Direct']
    assert vectortest(poscar.atoms[0],  [0.4094963891330107, 0.5, 0.5152348732651212])
    assert vectortest(poscar.atoms[2],  [0.17007786581988016, 0.5, 0.6638264732798937])
    assert vectortest(poscar.atoms[4],  [0.5905036108669895, 0.5, 0.48476512673487854])
    assert vectortest(poscar.atoms[6],  [0.8299221341801198, 0.5, 0.336173526720106])
    assert vectortest(poscar.atoms[8],  [0.38715052866080457, 0.5, 0.7846147360867455])
    assert vectortest(poscar.atoms[10], [0.11284947133919565, 0.5, 0.7846147360867455])
    assert vectortest(poscar.atoms[12], [0.25, 0.5, 0.36276037406667555])