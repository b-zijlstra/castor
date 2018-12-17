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
        self.outcar    = "OUTCAR"
        self.numbers   = None
        self.element   = None
        self.mass      = 2.0
        self.hessian   = "read"
        self.printmode = "normal"
        self.skip      = None
        self.first     = None
        self.temp      = [800]
        self.nummap    = None
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
                        print("Mass must be a number!")
                        sys.exit()
                readmode = None
                continue
            if(readmode=="nummap"):
                self.nummap = arg_in[i]
                readmode = None
                continue
            if(readmode=="skip"):
                self.skip = arg_in[i]
                readmode = None
                continue
            if(readmode=="first"):
                while True:
                    try:
                        self.first = int(arg_in[i])
                        break
                    except ValueError:
                        print("-f " + arg_in[i] + " must be an integer!")
                        sys.exit()
                readmode = None
                continue
            if(readmode=="temp"):
                self.temp = []
                templist = arg_in[i].split(',')
                for temperature in templist:
                    try:
                        self.temp.append(float(temperature))
                    except ValueError:
                        print(temperature + " must be a float!")
                        sys.exit()
                if(self.temp == 0):
                    print("self.temp is empty!")
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
            elif(sys.argv[i] == "-e" or sys.argv[i] == "--element"):
                readmode = "element"
                continue
            elif(sys.argv[i] == "-m" or sys.argv[i] == "--mass"):
                readmode = "mass"
                continue
            elif(sys.argv[i] == "--nummap"):
                readmode = "nummap"
                continue
            elif(sys.argv[i] == "-s" or sys.argv[i] == "--skip"):
                readmode = "skip"
                continue
            elif(sys.argv[i] == "-f" or sys.argv[i] == "--first"):
                readmode = "first"
                continue
            elif(sys.argv[i] == "-t" or sys.argv[i] == "--temp"):
                readmode = "temp"
                continue
            elif(sys.argv[i] == "-c" or sys.argv[i] == "--calc"):
                self.hessian = "calc"
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
                print("                   |    :|\n                   |     |\n                   |    .|\n               ____|    .|\n             .' .  ).   ,'\n           .' c   '7 ) (       nom-nom-nom\n       _.-\"       |.'   `.\n     .'           \"8E   :|\n     |          _}\"\"    :|\n     |         (   |     |\n    .'         )   |    :|\n/.beVER_.---.__8E  |    .|\n`BEver\"\"       \"\"  `-...-'")
                sys.exit()
            elif(sys.argv[i] == "--bahnhof"):
                print("  ____        _           _            __   _ \n |  _ \      | |         | |          / _| | |\n | |_) | __ _| |__  _ __ | |__   ___ | |_  | |\n |  _ < / _` | '_ \| '_ \| '_ \ / _ \|  _| | |\n | |_) | (_| | | | | | | | | | | (_) | |   |_|\n |____/ \__,_|_| |_|_| |_|_| |_|\___/|_|   (_)")
                sys.exit()
            else:
                print("Unexpected argument: " + sys.argv[i])
                self.help()
                sys.exit()
    def help(self):
        print("Use: getfreq.py <options>")
        print("Options:")
        print("-o or --outcar <outcar name>  | OUTCAR name to read. Example: $getfreq.py -o out (Default = OUTCAR)")
        print("-e or --element <element>     | Element type to change mass. Example: $getfreq.py -e H (Default = None)")
        print("-n or --numbers <numbers>     | Atom numbers of type. Example: $getfreq.py -n 34,36,38-40 (Default = all)")
        print("-m or --mass <mass>           | Mass to set for selected atoms. Example: $getfreq.py -m 13.0 (Default = 2.0)")
        print("--nummap <num=mass>           | Set atom number to certain mass. Example: $getfreq.py --nummap 34=13.0,35=2.0 (Default = None)")
        print("-s or --skip <numbers>        | Removes matrix elements for <numbers>. Example: $getfreq.py -s 38-40")
        print("-f or --first <number>        | Only print the first <number> of frequencies (Default = all)")
        print("-t or --temp <temperatures>   | Set temperatures for calculating partition functions. Example: $getfreq.py -t 493.15,513.15 (Default = 800 K)")
        print("-c or --calc                  | Calculate Hessian from forces. (Default = read Hessian from OUTCAR)")
        print("-l or --less                  | Set printmode to 'less' (Default = normal)")
        print("-ll or --least                | Set printmode to 'least' (Default = normal)")
        print("-a or --all                   | Set printmode to 'all' to also print displacements (Default = normal)")
        print("")
        print("-h or --help                  | displays this help message")

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    hessian = Hessian(arguments.temp)
    hessian.read(arguments.outcar,arguments.hessian)
    hessian.matrix.mass2diag()
    hessian.matrix.diag2freq(hessian.massmap)
    if(hessian.idipol > 0):
        hessian.getDipols()
        # hessian.writeDipols()
    if(arguments.nummap == None):
        hessian.mapMass(arguments.element, arguments.mass, arguments.numbers) # from numberlist, change all element masses to mass
    else:
        hessian.mapMass_nummap(arguments.nummap) #change mass of specific atom numbers
    hessian.setSkip(arguments.skip) # from skiplist, change all matrix elements to zero
    if(hessian.changes == True):
        hessian.addmatrix()
        hessian.newmatrices[0].sym2mass(hessian.massmap,hessian.skipset)
        hessian.newmatrices[0].mass2diag()
        hessian.newmatrices[0].diag2freq(hessian.massmap)
        if(hessian.idipol > 0):
            hessian.getNewIntens()
    hessian.write(arguments.printmode,arguments.first)

#EXECUTION
main(sys.argv)