#! /usr/bin/env python

# 
# getdis.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2017 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys

#MY PATHS

#MY CLASSES
from classes.class_outcar import Outcar

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.outcar    = "OUTCAR"
        self.numbers   = []
        self.printmode = "all"
    def setup(self, arg_in):
        readmode = None
        for i in range(1,len(arg_in)):
            if(readmode=="outcar"):
                self.outcar = arg_in[i]
                readmode = None
                continue
            if(readmode=="numbers"):
                numberstring = arg_in[i]
                numberpair = numberstring.split('-')
                if(len(numberpair) != 2):
                    print "StringError - invalid numberpair:" + numberstring
                    self.help()
                    sys.exit()
                try:
                    self.numbers = [int(numberpair[0]),int(numberpair[1])]
                except ValueError:
                    print "ValueError - invalid numberpair:" + numberstring
                    self.help()
                    sys.exit()
                readmode = None
                continue
            if(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-o" or sys.argv[i] == "--outcar"):
                readmode = "outcar"
                continue
            elif(sys.argv[i] == "-n" or sys.argv[i] == "--numbers"):
                readmode = "numbers"
                continue
                self.hessian = "calc"
                continue
            elif(sys.argv[i] == "-l" or sys.argv[i] == "--less"):
                self.printmode = "less"
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
        if(len(self.numbers) != 2):
            print "No atom pair selected. Did you specify '-n <numbers>'?"
            self.help()
            sys.exit()
    def help(self):
        print "Use: getdis.py <options>"
        print "Options:"
        print "-n or --numbers <numbers>     | Atom pair to get distance from. Example: $getdis.py -n 34-36 (No Default Value)"
        print "-o or --outcar <outcar name>  | OUTCAR name to read. Example: $getdis.py -o out (Default = OUTCAR)"
        print "-l or --less                  | Set printmode to 'less' to print only the last distance (Default = all)"
        print "-a or --all                   | Set printmode to 'all' to print the distance for all frames (Default = all)"
        print ""
        print "-h or --help                  | displays this help message"

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    outcar = Outcar()
    outcar.read(arguments.outcar)
    outcar.getDistance(arguments.numbers[0],arguments.numbers[1],arguments.printmode)

#EXECUTION
main(sys.argv)