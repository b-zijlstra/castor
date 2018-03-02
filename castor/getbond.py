#! /usr/bin/env python

# 
# getbond.py
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
from classes.class_chargemol import Chargemol

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.input    = "VASP_DDEC_analysis.output"
        self.atoms1   = []
        self.atoms2   = []
        self.printmode = "normal"
    def setup(self, arg_in):
        readmode = None
        for i in range(1,len(arg_in)):
            if(readmode=="input"):
                self.input = arg_in[i]
                readmode = None
                continue
            if(readmode=="atoms1"):
                numberstring = arg_in[i]
                numberstring = numberstring.split(',')
                try:
                    for number in numberstring:
                        self.atoms1.append(int(number))
                except ValueError:
                    print "ValueError - invalid atoms1:" + numberstring
                    self.help()
                    sys.exit()
                readmode = None
                continue
            if(readmode=="atoms2"):
                numberstring = arg_in[i]
                numberstring = numberstring.split(',')
                try:
                    for number in numberstring:
                        self.atoms2.append(int(number))
                except ValueError:
                    print "ValueError - invalid atoms1:" + numberstring
                    self.help()
                    sys.exit()
                readmode = None
                continue
            if(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-i" or sys.argv[i] == "--input"):
                readmode = "input"
                continue
            elif(sys.argv[i] == "-a1" or sys.argv[i] == "--atoms1"):
                readmode = "atoms1"
                continue
            elif(sys.argv[i] == "-a2" or sys.argv[i] == "--atoms2"):
                readmode = "atoms2"
                continue
            elif(sys.argv[i] == "-l" or sys.argv[i] == "--less"):
                self.printmode = "less"
                continue
            elif(sys.argv[i] == "-ll" or sys.argv[i] == "--least"):
                self.printmode = "least"
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
        print "Use: getbond.py <options>"
        print "Options:"
        print "-i or --input <input path>    | Input to read. Example: $getbond.py -i chargemol (Default = VASP_DDEC_analysis.output)"
        print "-a1 or --atoms1 <atoms1>      | Atoms to display bond orders of. Example: $getbond.py -a1 48,49 (Default = None)"
        print "-a2 or --atoms2 <atoms2>      | Specific atoms to connect to atoms of interest. Example: $getbond.py -a2 50,51 (Default = None)"
        print "-l or --less                  | Set printmode to 'less' (Default = normal)"
        print "-ll or --least                | Set printmode to 'least' (Default = normal)"
        print "-a or --all                   | Set printmode to 'all' to also print displacements (Default = normal)"
        print ""
        print "-h or --help                  | displays this help message"

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    chargemol = Chargemol()
    chargemol.read(arguments.input)
    chargemol.getbond(arguments.atoms1,arguments.atoms2,arguments.printmode)

#EXECUTION
main(sys.argv)