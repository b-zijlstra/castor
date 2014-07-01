#! /usr/bin/env python

# 
# test_getfreq.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 

# #IMPORTS
import sys, os, StringIO

#MY PATHS

#MY CLASSES
from castor.classes.class_hessian import Hessian

#DEFINES
def gethessian():
    hessian = Hessian()
    hessian.read(os.path.dirname(os.path.realpath(__file__)) + "/source/getfreq/OUTCAR")
    hessian.mapMass("H", 2.0, None) # from numberlist, change all element masses to mass
    hessian.newMassMatrix()
    hessian.diagonalize()
    hessian.calcFreqs()
    return hessian

def getfreq(hessian_in):
    output = StringIO.StringIO()
    sys.stdout = output
    hessian_in.write("all")
    sys.stdout = sys.__stdout__
    return output.getvalue()

def outputtest(output_in, file_in):
    # print output_in
    with open(file_in, 'r') as inputfile:
        lines = inputfile.readlines()
    output = output_in.split('\n')
    return all(x.rstrip() == y for x,y in zip(lines, output[4:]))

def vectortest(vec, check):
    print [vec.x, vec.y, vec.z]
    return all(round(x-y, 12) == 0 for x,y in zip([vec.x, vec.y, vec.z], check))

def test_hessian():
    hessian = gethessian()
    freq = getfreq(hessian)
    assert outputtest(freq, os.path.dirname(os.path.realpath(__file__)) + "/source/getfreq/freq")