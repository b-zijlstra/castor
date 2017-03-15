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
    count = 0
    summation = 1.0
    kb = 8.6173423e-2
    temp = 800.0
    for argument in arg_in:
        if count > 0:
            try:
                value = float(argument)
            except ValueError:
                count += 1
                continue
            if value != 0:
                summation *= 1.0 / (1.0 - exp(-value / (kb * temp)))
        count += 1
    print '{0:.{width}E}'.format(summation, width=6)

#EXECUTION
main(sys.argv)