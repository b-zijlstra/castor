#! /usr/bin/env python

# 
# calcnu.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2017 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
from math import exp

#MY PATHS

#MY CLASSES

#CLASSES

#DEFINES
def main(arg_in):
    summation = 1.0
    kb = 8.6173423e-2
    temp = 800
    if(len(arg_in)<=2):
        print("INCOMPLETE CALCNU.PY ARGUMENTS")
    else:
        try:
            temp = float(arg_in[1])
        except ValueError:
            temp = 800
            print("PROBLEM WITH CALCNU.PY TEMP ARGUMENT - RESETTING TO 800 K")
        for argument in arg_in[2:]:
            try:
                value = float(argument)
            except ValueError:
                continue
            if value != 0:
                summation *= 1.0 / (1.0 - exp(-value / (kb * temp)))
        print('{0:.{width}E}'.format(summation, width=6))

#EXECUTION
main(sys.argv)