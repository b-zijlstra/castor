#! /usr/bin/env python

# 
# xyz2pos.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2016 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys

#MY PATHS

#MY CLASSES
from classes.class_trajectory import Trajectory
from classes.class_poscar import Poscar

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.xyz = "out.xyz"
        self.posfile = "UNICAR"
        self.posmode = "Direct"
        self.elnames = False
    def setup(self, arg_in):
        self.readmode = None
        for i in range(1,len(arg_in)):
            if(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-c" or sys.argv[i] == "--car" or sys.argv[i] == "--cartesian"):
                self.posmode = "Cartesian"
                continue
            elif(sys.argv[i] == "-d" or sys.argv[i] == "--dir" or sys.argv[i] == "--direct"):
                self.posmode = "Direct"
                continue
            elif(sys.argv[i] == "-e" or sys.argv[i] == "--ele" or sys.argv[i] == "--elements"):
                self.elnames = True
                continue
            elif(sys.argv[i] == "-x" or sys.argv[i] == "--xyz"):
                self.readmode = "xyz"
                continue
            elif(sys.argv[i] == "-p" or sys.argv[i] == "--pos" or sys.argv[i] == "--poscar"):
                self.readmode = "poscar"
                continue
            elif(sys.argv[i] == "--bever"):
                print "                   |    :|\n                   |     |\n                   |    .|\n               ____|    .|\n             .' .  ).   ,'\n           .' c   '7 ) (       nom-nom-nom\n       _.-\"       |.'   `.\n     .'           \"8E   :|\n     |          _}\"\"    :|\n     |         (   |     |\n    .'         )   |    :|\n/.beVER_.---.__8E  |    .|\n`BEver\"\"       \"\"  `-...-'"
                sys.exit()
            elif(sys.argv[i] == "--bahnhof"):
                print "  ____        _           _            __   _ \n |  _ \      | |         | |          / _| | |\n | |_) | __ _| |__  _ __ | |__   ___ | |_  | |\n |  _ < / _` | '_ \| '_ \| '_ \ / _ \|  _| | |\n | |_) | (_| | | | | | | | | | | (_) | |   |_|\n |____/ \__,_|_| |_|_| |_|_| |_|\___/|_|   (_)"
                sys.exit()
            elif(self.readmode=="xyz"):
                self.xyz = arg_in[i]
                self.readmode = None
                continue
            elif(self.readmode=="poscar"):
                self.posfile = arg_in[i]
                self.readmode = None
                continue
            else:
                print "Unexpected argument: " + sys.argv[i]
                self.help()
                sys.exit()
    def help(self):
        print "Use: xyz2pos.py <options>"
        print "Options:"
        print "-x <xyz name>                     | Which xyz file to read. Example: $xyz2pos.py -x tot.xyz (Default = out.xyz)"
        print "-p <poscar file>                  | Which example poscar to use.  Example: $xyz2pos.py -p unitcell (Default = UNICAR)"
        print "-c or --car or --cartesian        | Fill POSCAR with cartesian coordinates. Example: $xyz2pos.py -c (Default = direct)"
        print "-d or --dir or --direct           | Fill POSCAR with direct coordinates. Example: $xyz2pos.py -d (Default = direct)"
        print "-e or --ele or --elements         | Include names of elements in POSCAR. Example: $xyz2pos.py -e (Default = False)"
        print ""
        print "-h or --help                      | displays this help message."

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    poscar = Poscar()
    poscar.read(arguments.posfile)
    trajectory = Trajectory()
    trajectory.readXYZ(arguments.xyz)
    if(arguments.elnames == True):
        poscar.elnames = trajectory.getLastElnames()
    poscar.elnr = trajectory.getLastElnr()
    poscar.extra = [arguments.posmode]
    poscar.atoms = trajectory.getLastCoord()
    if(arguments.posmode == "Direct"):
        poscar.unitcell.invertMatrix()
        poscar.cartesian2direct()
    poscar.write(rewrite=True)


#EXECUTION
main(sys.argv)