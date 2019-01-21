#! /usr/bin/env python3

# 
# class_rdf.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math

#MY PATHS

#MY CLASSES

class Rdf:
    """Defines a radial distribution function"""
    def __init__(self, trajectory_in, CNrange_in = 2.8):
        self.trajectory = trajectory_in
        self.bins = []
        self.distribution = []
        self.coordination = []
        self.CNrange = CNrange_in

    def genBins(self, binnumber_in = 100, binrange_in = 20.0):
        self.binnumber = binnumber_in
        self.binrange = float(binrange_in)
        self.delta = self.binrange/(self.binnumber); #binsize
        for i in range(0,self.binnumber):
            self.bins.append((i + 0.5) * self.delta)

        self.Vtot = self.trajectory.unitcell.getVolume()
        self.Vslice = 4.0/3.0 * math.pi * self.delta**3

    def sampleTrajectory(self, frame_start_in = None, frame_end_in = None):
        if(frame_start_in == None):
            frame_start = 0
        else:
            frame_start = int(frame_start_in) - 1

        if(frame_end_in == None):
            frame_end = len(self.trajectory.frames)
        else:
            frame_end = int(frame_end_in)  

        for i in range(frame_start,frame_end):
            self.sampleFrame(self.trajectory.frames[i])

    def sampleFrame(self,frame):
        dist = [];
        CN = [];

        for i in range(0,self.binnumber):
            dist.append(0.0)

        atoms = frame.atoms
        self.Npairs_inv = 1.0 / ( len(atoms) * (len(atoms)-1) )

        for i in range(0,len(atoms)):
            CN.append(0.0)

        for i in range(0,len(atoms)):
            for j in range(0,len(atoms)):
                if(j>i):
                    diffvec = atoms[j].vec - atoms[i].vec
                    if(self.trajectory.coordinates == "cartesian"):
                        diffvec = self.trajectory.unitcell.getModularCartesian(diffvec)
                        r2 = diffvec.length2()
                    elif(self.trajectory.coordinates == "direct"):
                        diffvec = self.trajectory.unitcell.getModularDirect(diffvec)
                        diffvec = self.trajectory.unitcell.direct2cartesian(diffvec)
                        r2 = diffvec.length2()
                    else:
                        print("Error: Unknown coordinate method")
                        sys.exit()

                    if(r2<self.CNrange**2):
                        CN[i] += 1
                        CN[j] += 1
                    if(r2<self.binrange**2):
                        dist[int(math.sqrt(r2)/self.delta)]+=2;

        self.coordination.append(CN)

        self.rdf = []
        for i in range(0,self.binnumber):
            Vi = self.Vslice * ((i+1)**3 - i**3)
            Vrel_inv = self.Vtot / Vi
            Drel = dist[i] * self.Npairs_inv # distribution /number of pairs
            self.rdf.append(Drel * Vrel_inv)
        self.distribution.append(self.rdf)

    def write(self, printmode_in = "all", rdfnum_in = None, decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces = spaces_in
        self.printmode = printmode_in

        if(self.printmode != "less"):
            print("RADIAL DISTRIBUTION")
        
        for i in range(0,len(self.bins)):
            string = ""

            x = '{0:.{width}f}'.format(self.bins[i], width=self.decimals)

            if(rdfnum_in==None):
                y = 0
                for j in range(0,len(self.distribution)):
                    y += self.distribution[j][i]
                y /= len(self.distribution)
                y = '{0:.{width}f}'.format(y, width=self.decimals)
            else:
                y = '{0:.{width}f}'.format(self.distribution[rdfnum_in][i], width=self.decimals)
            
            string += '{0:>{width}}'.format(x, width=self.decimals+self.spaces+3)
            string += '{0:>{width}}'.format(y, width=self.decimals+self.spaces+3)
            print(string)

        if(self.printmode != "less"):
            print("")
            print("COORDINATION NUMBERS")
            for i in range(0,17):
                string = '{0:>{width}}'.format(i, width=3)
                if(rdfnum_in==None):
                    count = 0.0
                    for j in range(0,len(self.coordination)):
                        for k in range(0,len(self.coordination[j])):
                            if(self.coordination[j][k] == i):
                                count += 1
                    count /= len(self.coordination)
                    string += '{0:>{width}}'.format(count, width=self.spaces+5)
                else:
                    count = 0
                    for k in range(0,len(self.coordination[rdfnum_in])):
                        if(self.coordination[rdfnum_in][k] == i):
                            count += 1
                    string += '{0:>{width}}'.format(count, width=self.spaces+5)
                print(string)


    def printTransient(self, printmode_in = "all", decimals_in = 6, spaces_in = 3):
        self.decimals = decimals_in
        self.spaces = spaces_in
        self.printmode = printmode_in

        if(self.printmode != "less"):
            print("RADIAL DISTRIBUTION")

        for i in range(0,len(self.bins)):
            string = ""

            x = '{0:.{width}f}'.format(self.bins[i], width=self.decimals)
            string += '{0:>{width}}'.format(x, width=self.decimals+self.spaces+3)

            for j in range(0,len(self.distribution)):
                y = self.distribution[j][i]
                y = '{0:.{width}f}'.format(y, width=self.decimals)
                string += '{0:>{width}}'.format(y, width=self.decimals+self.spaces+3)
            print(string)


        if(self.printmode != "less"):
            print("")
            print("COORDINATION NUMBERS")
            for i in range(0,17):
                string = '{0:>{width}}'.format(i, width=3)
                for j in range(0,len(self.coordination)):
                    count = 0
                    for k in range(0,len(self.coordination[j])):
                        if(self.coordination[j][k] == i):
                            count += 1
                    string += '{0:>{width}}'.format(count, width=self.spaces+5)
                print(string)