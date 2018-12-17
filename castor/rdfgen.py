#! /usr/bin/env python

# 
# rdfgen.py
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
from classes.class_trajectory import Trajectory
from classes.class_rdf import Rdf
from classes.class_rdfanalyser import RdfAnalyser

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.xyz = []
        self.bins = 100
        self.binrange = 20
        self.multimode = None
        self.printmode = "all"
        self.periodic = [False, False, False]
        self.CNrange = 2.8
    def setup(self, arg_in):
        self.readmode = None
        for i in range(1,len(arg_in)):
            if(self.readmode=="bins"):
                self.bins = int(arg_in[i])
                self.readmode = None
                continue
            elif(self.readmode=="binrange"):
                self.binrange = float(arg_in[i])
                self.readmode = None
                continue
            elif(self.readmode=="periodic"):
                s = arg_in[i]
                if(s.find('x')!=-1 or s.find('X')!=-1):
                    self.periodic[0] = True
                if(s.find('y')!=-1 or s.find('Y')!=-1):
                    self.periodic[1] = True
                if(s.find('z')!=-1 or s.find('Z')!=-1):
                    self.periodic[2] = True
                self.readmode = None
                continue
            elif(self.readmode=="CNrange"):
                self.CNrange = float(arg_in[i])
                self.readmode = None
                continue
            elif(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-x" or sys.argv[i] == "--xyz"):
                self.readmode = "xyz"
                continue
            elif(sys.argv[i] == "-b" or sys.argv[i] == "--mb" or sys.argv[i] == "--bins"):
                self.readmode = "bins"
                continue
            elif(sys.argv[i] == "-r" or sys.argv[i] == "--br" or sys.argv[i] == "--range" or sys.argv[i] == "--binrange"):
                self.readmode = "binrange"
                continue
            elif(sys.argv[i] == "-p" or sys.argv[i] == "--periodic"):
                self.readmode = "periodic"
                continue
            elif(sys.argv[i] == "--CNrange" or sys.argv[i] == "--crange"):
                self.readmode = "CNrange"
                continue
            elif(sys.argv[i] == "-l" or sys.argv[i] == "--less"):
                self.printmode = "less"
                self.readmode = None
                continue
            elif(sys.argv[i] == "--overlay"):
                self.multimode = "overlay"
                self.readmode = None
                continue
            elif(sys.argv[i] == "--trans" or sys.argv[i] == "--transient"):
                self.multimode = "transient"
                self.readmode = None
                continue
            elif(sys.argv[i] == "--average" or sys.argv[i] == "--avg"):
                self.multimode = "average"
                self.readmode = None
                continue
            elif(sys.argv[i] == "--bever"):
                print("                   |    :|\n                   |     |\n                   |    .|\n               ____|    .|\n             .' .  ).   ,'\n           .' c   '7 ) (       nom-nom-nom\n       _.-\"       |.'   `.\n     .'           \"8E   :|\n     |          _}\"\"    :|\n     |         (   |     |\n    .'         )   |    :|\n/.beVER_.---.__8E  |    .|\n`BEver\"\"       \"\"  `-...-'")
                sys.exit()
            elif(sys.argv[i] == "--bahnhof"):
                print("  ____        _           _            __   _ \n |  _ \      | |         | |          / _| | |\n | |_) | __ _| |__  _ __ | |__   ___ | |_  | |\n |  _ < / _` | '_ \| '_ \| '_ \ / _ \|  _| | |\n | |_) | (_| | | | | | | | | | | (_) | |   |_|\n |____/ \__,_|_| |_|_| |_|_| |_|\___/|_|   (_)")
                sys.exit()
            elif(self.readmode=="xyz"):
                self.xyz.append(arg_in[i])
                continue
            else:
                print("Unexpected argument: " + sys.argv[i])
                self.help()
                sys.exit()
    def help(self):
        print("Use: rdfgen.py <options>")
        print("Options:")
        print("-x <xyz name> <xyz_2 name>       | Example: $rdfgen.py -o fcc.xyz (Default = out.xyz)")
        print("-b or --mb or --bins               | Sets the maximum number of bins for the rdf. Example: $rdfgen.py -b 50 (Default = 100)")
        print("-r or --br or --range or --binrange | Sets the bin range for the rdf in angstrom. Example: $rdfgen.py -r 10 (Default = 20)")
        print("-l or --less                      | Only the most important output is printed.")
        print("--overlay                         | When multiple trajectories are given, use same bins and print an rdf column for each averaged trajectory.")
        print("--avg or --average                 | When multiple trajectories are given, average all frames and print one rdf column.")
        print("--trans or --transient             | For all trajectories use overlay for every frame.")
        print("")
        print("-h or --help                      | displays this help message.")

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    rdflist = []
    for xyz in arguments.xyz:
        trajectory = Trajectory()
        trajectory.readXYZ(xyz)
        trajectory.setBoundaries(*arguments.periodic)
        rdf = Rdf(trajectory, arguments.CNrange)
        rdf.genBins(arguments.bins, arguments.binrange)
        rdf.sampleTrajectory()
        rdflist.append(rdf)
    if(arguments.multimode == "average"):
        analyse = RdfAnalyser(rdflist)
        analyse.checkSimilarity()
        analyse.printAverage(arguments.printmode)
    elif(arguments.multimode == "overlay"):
        analyse = RdfAnalyser(rdflist)
        analyse.checkSimilarity()
        analyse.printOverlay(arguments.printmode)
    elif(arguments.multimode == "transient"):
        for rdf in rdflist:
            rdf.printTransient(arguments.printmode)
            print("")
    else:
        for rdf in rdflist:
            rdf.write(arguments.printmode)
            print("")



#EXECUTION
main(sys.argv)