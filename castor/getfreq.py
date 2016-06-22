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
        self.numbers = None
        self.element = "H"
        self.mass = 2.0
        self.printmode = "all"
        self.skip = False
    def setup(self, arg_in):
        readmode = None
        for i in range(1,len(arg_in)):
            if(readmode=="outcar"):
                self.outcar = arg_in[i]
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
            if(sys.argv[i] == "-s" or sys.argv[i] == "--skip"):
                self.skip = True
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
            elif(sys.argv[i] == "-e" or sys.argv[i] == "--element"):
                readmode = "element"
                continue
            elif(sys.argv[i] == "-m" or sys.argv[i] == "--mass"):
                readmode = "mass"
                continue
            elif(sys.argv[i] == "-l" or sys.argv[i] == "--less"):
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
        print "-o or --outcar <outcar name>  | Example: $getfreq.py -o out (Default = OUTCAR)"
        print "-n or --numbers <numbers>     | Example: $getfreq.py -n 34,36,38-40 (Default = all)"
        print "-e or --element <element>     | Example: $getfreq.py -e C (Default = H)"
        print "-m or --mass <mass>           | Example: $getfreq.py -m 13.0 (Default = 2.0)"
        print "-s or --skip                  | Removes matrix elements for -n. Example: $getfreq.py -s -n 38-40"
        print "-l or --less                  | Set printmode to 'less' (Default = all)"
        print ""
        print "-h or --help                  | displays this help message"

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    hessian = Hessian()
    hessian.read(arguments.outcar)
    hessian.matrix.mass2diag()
    hessian.matrix.diag2freq(hessian.massmap)
    if(hessian.idipol > 0):
        hessian.getDipols()
        hessian.writeDipols()
    if(arguments.skip == False):
        hessian.mapMass(arguments.element, arguments.mass, arguments.numbers) # from numberlist, change all element masses to mass
        hessian.setSkip(None)
    elif(arguments.skip == True):
        hessian.mapMass(None, None, None)
        hessian.setSkip(arguments.numbers) # from numberlist, change all matrix elements to zero
    if(hessian.changes == True):
        hessian.addmatrix()
        hessian.newmatrices[0].sym2mass(hessian.massmap,hessian.skipset)
        hessian.newmatrices[0].mass2diag()
        hessian.newmatrices[0].diag2freq(hessian.massmap)
        if(hessian.idipol > 0):
            hessian.getNewIntens()
    hessian.write(arguments.printmode)

#EXECUTION
main(sys.argv)