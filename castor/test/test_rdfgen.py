#! /usr/bin/env python

# 
# test_rdf.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 

# #IMPORTS
import sys, os, StringIO

#MY PATHS

#MY CLASSES
from castor.classes.class_trajectory import Trajectory
from castor.classes.class_rdf import Rdf
from castor.classes.class_rdfanalyser import RdfAnalyser

#DEFINES
def getanalyse():
    rdflist = []
    for xyz in [os.path.dirname(os.path.realpath(__file__)) + "/source/rdfgen/test1.xyz", os.path.dirname(os.path.realpath(__file__)) + "/source/rdfgen/test2.xyz"]:
        trajectory = Trajectory()
        trajectory.readXYZ(xyz)
        trajectory.setBoundaries(True, True, True)
        rdf = Rdf(trajectory, 2.8)
        rdf.genBins(10, 10)
        rdf.sampleTrajectory()
        rdflist.append(rdf)
    analyse = RdfAnalyser(rdflist)
    analyse.checkSimilarity()
    return analyse

def getAverage(analyse_in):
    output = StringIO.StringIO()
    sys.stdout = output
    analyse_in.printAverage("all")
    sys.stdout = sys.__stdout__
    return output.getvalue()

def getOverlay(analyse_in):
    output = StringIO.StringIO()
    sys.stdout = output
    analyse_in.printOverlay("all")
    sys.stdout = sys.__stdout__
    return output.getvalue()

def outputtest(output_in, file_in):
    # print output_in
    with open(file_in, 'r') as inputfile:
        lines = inputfile.read()
    return output_in == lines

def test_analyse():
    analyse = getanalyse()
    average = getAverage(analyse)
    overlay = getOverlay(analyse)
    assert outputtest(average, os.path.dirname(os.path.realpath(__file__)) + "/source/rdfgen/average")
    assert outputtest(overlay, os.path.dirname(os.path.realpath(__file__)) + "/source/rdfgen/overlay")