#! /usr/bin/env python

# 
# unitcopy.py
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
        self.unifile = "POSCAR1"
        self.posfile = "POSCAR2"
    def setup(self, arg_in):
        self.readmode = None
        for i in range(1,len(arg_in)):
            if(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()

            elif(sys.argv[i] == "-u" or sys.argv[i] == "--uni" or sys.argv[i] == "--unitcell"):
                self.readmode = "unifile"
                continue
            elif(sys.argv[i] == "-p" or sys.argv[i] == "--pos" or sys.argv[i] == "--poscar"):
                self.readmode = "poscar"
                continue
            elif(sys.argv[i] == "--bever"):
                print("                   |    :|\n                   |     |\n                   |    .|\n               ____|    .|\n             .' .  ).   ,'\n           .' c   '7 ) (       nom-nom-nom\n       _.-\"       |.'   `.\n     .'           \"8E   :|\n     |          _}\"\"    :|\n     |         (   |     |\n    .'         )   |    :|\n/.beVER_.---.__8E  |    .|\n`BEver\"\"       \"\"  `-...-'")
                sys.exit()
            elif(sys.argv[i] == "--bahnhof"):
                print("  ____        _           _            __   _ \n |  _ \      | |         | |          / _| | |\n | |_) | __ _| |__  _ __ | |__   ___ | |_  | |\n |  _ < / _` | '_ \| '_ \| '_ \ / _ \|  _| | |\n | |_) | (_| | | | | | | | | | | (_) | |   |_|\n |____/ \__,_|_| |_|_| |_|_| |_|\___/|_|   (_)")
                sys.exit()
            elif(self.readmode=="unifile"):
                self.unifile = arg_in[i]
                self.readmode = None
                continue
            elif(self.readmode=="poscar"):
                self.posfile = arg_in[i]
                self.readmode = None
                continue
            else:
                print("Unexpected argument: " + sys.argv[i])
                self.help()
                sys.exit()
    def help(self):
        print("Use: unitcopy.py <options>")
        print("Options:")
        print("-u <unitcell file>                | Which poscar to use for unitcell. Example: $unitcopy.py -u unitcell (Default = POSCAR1)")
        print("-p <poscar file>                  | Which poscar to use for positions. Example: $unitcopy.py -p poscar (Default = POSCAR2)")
        print("")
        print("-h or --help                      | displays this help message.")

#DEFINES

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    poscar_unit = Poscar()
    poscar_pos = Poscar()
    poscar_unit.read(arguments.unifile)
    poscar_pos.read(arguments.posfile)
    poscar_pos.unitcell = poscar_unit.unitcell
    poscar_pos.write(rewrite=True)


#EXECUTION
main(sys.argv)