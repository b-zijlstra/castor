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
def gethessian(file_in):
    hessian = Hessian()
    hessian.read(file_in)
    hessian.mapMass("H", 2.0, None) # from numberlist, change all element masses to mass
    hessian.newMassMatrix()
    hessian.diagonalize()
    hessian.calcFreqs()
    return hessian

def getfreq(hessian_in, writemode_in = "all"):
    output = StringIO.StringIO()
    sys.stdout = output
    hessian_in.write(writemode_in)
    sys.stdout = sys.__stdout__
    return output.getvalue()

def outputtest(output_in, file_in, offset_in = 0):
    # print output_in
    with open(file_in, 'r') as inputfile:
        lines = inputfile.readlines()
    output = output_in.split('\n')
    for x,y in zip(lines, output[offset_in:]):
        assert x.rstrip() == y 
    

def vectortest(vec, check):
    print [vec.x, vec.y, vec.z]
    return all(round(x-y, 12) == 0 for x,y in zip([vec.x, vec.y, vec.z], check))

def test_hessian1():
    hessian = gethessian(os.path.dirname(os.path.realpath(__file__)) + "/source/getfreq/OUTCAR")
    freq = getfreq(hessian)
    outputtest(freq, os.path.dirname(os.path.realpath(__file__)) + "/source/getfreq/freq", 4)

def test_hessian2():
    hessian = gethessian(os.path.dirname(os.path.realpath(__file__)) + "/source/getfreq/OUTCAR2")
    freq = getfreq(hessian, "less")
    outputtest(freq, os.path.dirname(os.path.realpath(__file__)) + "/source/getfreq/freq2")