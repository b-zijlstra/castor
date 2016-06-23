#! /usr/bin/env python

# 
# intense.py
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
from classes.class_dipol import Dipols

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.outcar    = "OUTCAR"
        self.printmode = "intensity"
    def setup(self, arg_in):
        readmode = None
        for i in range(1,len(arg_in)):
            if(readmode=="outcar"):
                self.outcar = arg_in[i]
                readmode = None
                continue
            if(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-o" or sys.argv[i] == "--outcar"):
                readmode = "outcar"
                continue
            elif(sys.argv[i] == "-a" or sys.argv[i] == "--all"):
                self.printmode = "all"
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
        print "Use: intense.py <options>"
        print "Options:"
        print "-o or --outcar <outcar name>  | Example: $intense.py -o out (Default = OUTCAR)"
        print "-a or --all                   | Print additional details. Example: $intense.py -a (Default = intensity)"
        print ""
        print "-h or --help                  | displays this help message"

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    dipols = Dipols()
    dipols.read(arguments.outcar)
    dipols.getMatrix()
    dipols.write(printmode=arguments.printmode)

#EXECUTION
main(sys.argv)