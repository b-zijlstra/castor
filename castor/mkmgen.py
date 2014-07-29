#! /usr/bin/env python

# 
# mkmgen.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys,subprocess

#MY PATHS

#MY CLASSES
from classes.class_hessian import Hessian

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.parfile = "par.mkm"
        self.sourcefile = "source.mkm"
        self.run = False
        self.getOutput = False
    def setup(self, arg_in):
        readmode = None
        for i in range(1,len(arg_in)):
            if(readmode=="parfile"):
                self.parfile = arg_in[i]
                readmode = None
                continue
            if(readmode=="sourcefile"):
                self.sourcefile = arg_in[i]
                readmode = None
                continue
            if(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-p" or sys.argv[i] == "--par"):
                readmode = "parfile"
                continue
            elif(sys.argv[i] == "-s" or sys.argv[i] == "--source"):
                readmode = "sourcefile"
                continue
            elif(sys.argv[i] == "-r" or sys.argv[i] == "--run"):
                self.run = True
                continue
            elif(sys.argv[i] == "-o" or sys.argv[i] == "--out"):
                self.getOutput = True
                continue
            elif(sys.argv[i] == "--bever"):
                print "                   |    :|\n                   |     |\n                   |    .|\n               ____|    .|\n             .' .  ).   ,'\n           .' c   '7 ) (       nom-nom-nom\n       _.-\"       |.'   `.\n     .'           \"8E   :|\n     |          _}\"\"    :|\n     |         (   |     |\n    .'         )   |    :|\n/.beVER_.---.__8E  |    .|\n`BEver\"\"       \"\"  `-...-'"
                sys.exit()
            elif(sys.argv[i] == "--bahnhof"):
                print "  ____        _           _            __   _ \n |  _ \      | |         | |          / _| | |\n | |_) | __ _| |__  _ __ | |__   ___ | |_  | |\n |  _ < / _` | '_ \| '_ \| '_ \ / _ \|  _| | |\n | |_) | (_| | | | | | | | | | | (_) | |   |_|\n |____/ \__,_|_| |_|_| |_|_| |_|\___/|_|   (_)"
                sys.exit()
            else:
                print "Unexpected argument: " + sys.argv[i]
                self.help()
                sys.exit()
    def help(self):
        print "Use: mkmgen.py <options>"
        print "Options:"
        print "-p or --par <parfile name>       | Example: $mkmgen.py -p parfile (Default = par.mkm)"
        print "-s or --source <source name>     | Example: $mkmgen.py -p sourcefile (Default = source.mkm)"
        print "-r or --run                      | Automatically run all created input."
        print "-o or --out                      | Automatically get all created output."
        print ""
        print "-h or --help                     | displays this help message"

#DEFINES

def makemkm(parfile_in, sourcefile_in):
    with open(parfile_in, 'r') as par:
        lines = par.readlines()
        for line in lines[1:]:
            line = line.split()
            subprocess.call(["mkmmod.sh", str(line[0]), str(line[1]), str(line[2]), str(line[3]), str(line[4]), str(line[5]), str(line[6]), str(line[7]), sourcefile_in])

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    makemkm(arguments.parfile, arguments.sourcefile)
    if(arguments.run == True):
        subprocess.call(["mkmrun.sh"])
    if(arguments.getOutput == True):
        subprocess.call(["mkmout.sh"])

#EXECUTION
main(sys.argv)