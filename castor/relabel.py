#! /usr/bin/env python

# 
# relabel.py
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
def main(input1_in, input2_in, target = "metal"):
	poscar2_in = Poscar()
	poscar_out = Poscar()
	poscar2_in.read(input2_in)
	poscar_out.read(input1_in)
	poscar_out.relabel(poscar2_in, target)
	poscar_out.write()

#EXECUTION
if len(sys.argv)==3:
	arguments = sys.argv[1], sys.argv[2], "metal"
	main(*arguments)
elif len(sys.argv)==4:
	arguments = sys.argv[1], sys.argv[2], sys.argv[3]
	main(*arguments)
else:
	print("Unexpected amount of arguments.")