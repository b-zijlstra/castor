#! /usr/bin/env python

# 
# getfreq.py
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
from classes.class_hessian import Hessian

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.outcar = "OUTCAR"
        self.poscar = "POSCAR"
        self.numbers = None
        self.element = "H"
        self.mass = 2.0
        self.printmode = "all"
    def setup(self, arg_in):
        readmode = None
        for i in range(1,len(arg_in)):
            if(readmode=="outcar"):
                self.outcar = arg_in[i]
                readmode = None
                continue
            if(readmode=="poscar"):
                self.poscar = arg_in[i]
                readmode = None
                continue
            if(readmode=="numbers"):
                self.numbers = arg_in[i]
                readmode = None
                continue
            if(readmode=="element"):
                self.element = arg_in[i]
                readmode = None
                continue
            if(readmode=="mass"):
                while True:
                    try:
                        self.mass = float(arg_in[i])
                        if(self.mass <= 0):
                            self.mass = None
                        break
                    except ValueError:
                        print "Mass must be a number!"
                        sys.exit()
                readmode = None
                continue

            if(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-o"):
                readmode = "outcar"
                continue
            elif(sys.argv[i] == "-p"):
                readmode = "poscar"
                continue
            elif(sys.argv[i] == "-n"):
                readmode = "numbers"
                continue
            elif(sys.argv[i] == "-e"):
                readmode = "element"
                continue
            elif(sys.argv[i] == "-m"):
                readmode = "mass"
                continue
            elif(sys.argv[i] == "--less"):
                self.printmode = "less"
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
    def getOutcar(self):
        return self.outcar
    def getPoscar(self):
        return self.poscar
    def getNumbers(self):
        return self.numbers
    def getElement(self):
        return self.element
    def getMass(self):
        return self.mass
    def getPrintmode(self):
        return self.printmode
    def help(self):
        print "Use: getfreq.py <options>"
        print "Options:"
        print "-o <outcar name>      | Example: $getfreq.py -o out (Default = OUTCAR)"
        print "-p <poscar name>      | Example: $getfreq.py -p pos (Default = POSCAR)"
        print "-n <numbers>          | Example: $getfreq.py -n 34,36,38-40 (Default = all)"
        print "-e <element>          | Example: $getfreq.py -e C (Default = H)"
        print "-m <mass>             | Example: $getfreq.py -m 13.0 (Default = 2.0)"
        print ""
        print "-h or --help           | displays this help message"

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    hessian = Hessian()
    hessian.read(arguments.getOutcar(), arguments.getPoscar())
    # hessian.mapMass() # do not changes masses
    # hessian.mapMass("H", 2.0) # change hydrogen mass to 2.0 (deuterium)
    hessian.mapMass(arguments.getElement(), arguments.getMass(), arguments.getNumbers()) # from numberlist, change all element masses to mass
    # hessian.mapMass("C", 13.003) # change carbon mass to carbon-13
    hessian.newMassMatrix()
    hessian.diagonalize()
    hessian.calcFreqs()
    hessian.write(arguments.getPrintmode())


#EXECUTION
main(sys.argv)