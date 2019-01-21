#! /usr/bin/env python3

# 
# rotate.py
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
from classes.class_poscar import Poscar
from classes.class_rotation import Matrix, AxisAngle
from classes.class_vector import Vector

#DEFINES
def main(rotX1_in, rotY1_in, rotZ1_in, rotX2_in, rotY2_in, rotZ2_in, angle_in, input_in):
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
	poscar.write()

#EXECUTION
if len(sys.argv)==8:
	arguments = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), "POSCAR"
	main(*arguments)
elif len(sys.argv)==9:
	arguments = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
	main(*arguments)
else:
	print("Unexpected amount of arguments.")