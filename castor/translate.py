#! /usr/bin/env python

# 
# translate.py
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

#DEFINES
def main(transX_in, transY_in, sizeX_in, sizeY_in, input_in):
	poscar_in = Poscar()
	poscar_out = Poscar()
	poscar_in.read(input_in)
	poscar_out.read(input_in)
	poscar_out.translate(transX_in, transY_in, sizeX_in, sizeY_in)
	poscar_out.move2unitcell()
	poscar_out.relabel(poscar_in,"metal")
	poscar_out.write()

#EXECUTION
if len(sys.argv)==5:
	arguments = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), "POSCAR"
	main(*arguments)
elif len(sys.argv)==6:
	arguments = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
	main(*arguments)
else:
	print("Unexpected amount of arguments.")