#! /usr/bin/env python

# 
# class_rdfanalyser.py
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

class RdfAnalyser:
    """Defines a RdfAnalyser, which can hold multiple rdfs"""
    def __init__(self, rdfs_in):
        self.rdfs = rdfs_in
        self.sameRange = False  #if true, all rdfs have similar x-range
        self.sameProperties = False  #if true, all rdfs have similar x-range and properties
    def checkSimilarity(self): #check if properties are the same so we can average/overlay rdfs
        self.sameRange = True
        self.sameProperties = True
        for rdf in self.rdfs:
            self.reference = self.rdfs[0]
            if(rdf.bins != self.reference.bins):
                self.sameRange = False
                self.sameProperties = False
                break

            a = rdf.trajectory.unitcell
            b = self.reference.trajectory.unitcell
            list_id = ["vec_1.x", "vec_1.y", "vec_1.z", "vec_2.x", "vec_2.y", "vec_2.z", "vec_3.x", "vec_3.y", "vec_3.z", "periodic_1", "periodic_2", "periodic_3"]
            list_a = [a.vec_1.x, a.vec_1.y, a.vec_1.z, a.vec_2.x, a.vec_2.y, a.vec_2.z, a.vec_3.x, a.vec_3.y, a.vec_3.z, a.periodic_1, a.periodic_2, a.periodic_3]
            list_b = [b.vec_1.x, b.vec_1.y, b.vec_1.z, b.vec_2.x, b.vec_2.y, b.vec_2.z, b.vec_3.x, b.vec_3.y, b.vec_3.z, b.periodic_1, b.periodic_2, b.periodic_3]
            for i in range(0,len(list_a)):
                if(list_a[i] != list_b[i]):
                    self.sameProperties = False
                    # print("RDF properties not the same for " + list_id[i] + ": " + str(list_a[i]) + " != " + str(list_b[i]))

    def printAverage(self, printmode_in = "all", decimals_in = 6, spaces_in = 3):
        if(self.sameRange==True):
            self.decimals = decimals_in
            self.spaces = spaces_in
            self.printmode = printmode_in

            if(self.printmode != "less"):
                print("RADIAL DISTRIBUTION")
            
            for i in range(0,len(self.reference.bins)):
                string = ""

                x = '{0:.{width}f}'.format(self.reference.bins[i], width=self.decimals)

                y = 0
                for j in range(0,len(self.rdfs)):
                    for k in range(0,len(self.rdfs[j].distribution)):
                        y += self.rdfs[j].distribution[k][i]
                y /= (len(self.rdfs) * len(self.rdfs[0].distribution))
                y = '{0:.{width}f}'.format(y, width=self.decimals)

                string += '{0:>{width}}'.format(x, width=self.decimals+self.spaces+3)
                string += '{0:>{width}}'.format(y, width=self.decimals+self.spaces+3)
                print(string)

            if(self.printmode != "less"):
                print("")
                print("COORDINATION NUMBERS")
                for i in range(0,17):
                    string = '{0:>{width}}'.format(i, width=3)
                    count = 0.0
                    for j in range(0,len(self.rdfs)):
                        for k in range(0,len(self.rdfs[j].coordination)):
                            for l in range(0,len(self.rdfs[j].coordination[k])):
                                if(self.rdfs[j].coordination[k][l] == i):
                                    count += 1
                    count /= (len(self.rdfs) * len(self.rdfs[0].coordination))
                    string += '{0:>{width}}'.format(count, width=self.spaces+5)
                    print(string)
        else:
            print("Can not average the rdfs. X-range is not the same!")

    def printOverlay(self, printmode_in = "all", decimals_in = 6, spaces_in = 3):
        if(self.sameRange==True):
            self.decimals = decimals_in
            self.spaces = spaces_in
            self.printmode = printmode_in

            if(self.printmode != "less"):
                print("RADIAL DISTRIBUTION")
            
            for i in range(0,len(self.reference.bins)):
                string = ""

                x = '{0:.{width}f}'.format(self.reference.bins[i], width=self.decimals)
                string += '{0:>{width}}'.format(x, width=self.decimals+self.spaces+3)

                for j in range(0,len(self.rdfs)):
                    y = 0
                    for k in range(0,len(self.rdfs[j].distribution)):
                        y += self.rdfs[j].distribution[k][i]
                    y /= len(self.rdfs[j].distribution)
                    y = '{0:.{width}f}'.format(y, width=self.decimals)
                    string += '{0:>{width}}'.format(y, width=self.decimals+self.spaces+3)
                print(string)

            if(self.printmode != "less"):
                print("")
                print("COORDINATION NUMBERS")
                for i in range(0,17):
                    string = '{0:>{width}}'.format(i, width=3)
                    for j in range(0,len(self.rdfs)):
                        count = 0.0
                        for k in range(0,len(self.rdfs[j].coordination)):
                            for l in range(0,len(self.rdfs[j].coordination[k])):
                                if(self.rdfs[j].coordination[k][l] == i):
                                    count += 1
                        string += '{0:>{width}}'.format(count, width=self.spaces+5)
                    print(string)

        else:
            print("Can not overlay the rdfs. X-range is not the same!")