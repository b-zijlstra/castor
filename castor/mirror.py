#! /usr/bin/env python

# 
# mirror.py
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
from classes.class_poscar import Poscar

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.posfile = "POSCAR"
        self.mode = ["all", "top"]
        self.symmetry = ["mirror", "mirror"] #Use 'mirror','cell_center', 'atom_center'
        self.zplane = 0.5 #direct middle of the z-axis of the structure
        self.zrange = 0.01 #deviation tolerance in angstrom
        self.getOffset = False
        self.setOffset = None
    def setup(self, arg_in):
        self.readmode = None
        for i in range(1,len(arg_in)):
            if(self.readmode=="symmetry"):
                s = arg_in[i].split('-')
                if(s[0].find('mir')!=-1):
                    self.symmetry[0] = "mirror"
                    self.symmetry[1] = "mirror"
                elif(s[0].find('cell')!=-1):
                    self.symmetry[0] = "cell_center"
                    self.symmetry[1] = "cell_center"
                elif(s[0].find('atom')!=-1):
                    self.symmetry[0] = "atom_center"
                    self.symmetry[1] = "atom_center"
                elif(s[0].find('off')!=-1):
                    self.symmetry[0] = "offset"
                    self.symmetry[1] = "offset"
                if(len(s)==2):
                    if(s[1].find('mir')!=-1):
                        self.symmetry[1] = "mirror"
                    elif(s[1].find('cell')!=-1):
                        self.symmetry[1] = "cell_center"
                    elif(s[1].find('atom')!=-1):
                        self.symmetry[1] = "atom_center"
                    elif(s[1].find('off')!=-1):
                        self.symmetry[1] = "offset"
                self.readmode = None
                continue
            elif(self.readmode=="zplane"):
                self.zplane = float(arg_in[i])
                self.readmode = None
                continue
            elif(self.readmode=="zrange"):
                self.zrange = float(arg_in[i])
                self.readmode = None
                continue
            elif(self.readmode=="posfile"):
                self.posfile = arg_in[i]
                self.readmode = None
                continue
            elif(self.readmode=="offset"):
                s = arg_in[i].split(':')
                if(len(s)==3):
                    self.setOffset = [s[0], s[1], s[2]]
                else:
                    print "Offset input has wrong length: 3 != " + str(len(s))
                    self.help()
                    sys.exit()
                self.readmode = None
                continue
            elif(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-s" or sys.argv[i] == "--sym" or sys.argv[i] == "--symmetry"):
                self.readmode = "symmetry"
                continue
            elif(sys.argv[i] == "-z" or sys.argv[i] == "--zp" or sys.argv[i] == "--zplane"):
                self.readmode = "zplane"
                continue
            elif(sys.argv[i] == "-r" or sys.argv[i] == "--zr" or sys.argv[i] == "--zrange"):
                self.readmode = "zrange"
                continue
            elif(sys.argv[i] == "-p" or sys.argv[i] == "--pos" or sys.argv[i] == "--poscar"):
                self.readmode = "posfile"
                continue
            elif(sys.argv[i] == "--setOffset"):
                self.readmode = "offset"
                continue
            elif(sys.argv[i] == "-a" or sys.argv[i] == "--ads" or sys.argv[i] == "--adsorbate"):
                self.mode[0] = "ads"
                self.readmode = None
                continue
            elif(sys.argv[i] == "-b" or sys.argv[i] == "--bottom"):
                self.mode[1] = "bottom"
                self.readmode = None
                continue
            elif(sys.argv[i] == "-t" or sys.argv[i] == "--top"):
                self.mode[1] = "top"
                self.readmode = None
                continue
            elif(sys.argv[i] == "-g" or sys.argv[i] == "--getOffset"):
                self.getOffset = True
                continue
            elif(sys.argv[i] == "-bever"):
                print "                   |    :|\n                   |     |\n                   |    .|\n               ____|    .|\n             .' .  ).   ,'\n           .' c   '7 ) (       nom-nom-nom\n       _.-\"       |.'   `.\n     .'           \"8E   :|\n     |          _}\"\"    :|\n     |         (   |     |\n    .'         )   |    :|\n/.beVER_.---.__8E  |    .|\n`BEver\"\"       \"\"  `-...-'"
                sys.exit()
            elif(sys.argv[i] == "-bahnhof"):
                print "  ____        _           _            __   _ \n |  _ \      | |         | |          / _| | |\n | |_) | __ _| |__  _ __ | |__   ___ | |_  | |\n |  _ < / _` | '_ \| '_ \| '_ \ / _ \|  _| | |\n | |_) | (_| | | | | | | | | | | (_) | |   |_|\n |____/ \__,_|_| |_|_| |_|_| |_|\___/|_|   (_)"
                sys.exit()
            else:
                print "Unexpected argument: " + sys.argv[i]
                self.help()
                sys.exit()
    def help(self):
        print "Use: mirror.py <options>"
        print "Options:"
        print "-s or --sym or --symmetry <sym>          | Set symmetry mode. Example: $mirror.py -s 'cell-atom' (Default = mirror-mirror)"
        print "-z or --zp or --zplane <val>             | Direct middle of the z-axis of the structure. Example: $mirror.py -z 0.6 (Default = 0.5)"
        print "-r or --zr or --zrange <range>           | Z-deviation tolerance in angstrom. Example: $mirror.py -r 1.0 (Default = 0.01)"
        print "-p or --pos or --poscar <POSCAR name>    | Example: $mirror.py -p CONTCAR (Default = POSCAR)"
        print "-a or --ads or --adsorbate               | Mirror only the adsorbate (Default = all)"
        print "-b or --bottom                           | Mirror bottom atoms to top instead of top to bottom (Default = top)"
        print "--setOffset                              | Set the offset vector Example: $mirror.py --setOffset 0.1:0.2:0.3 (Default = None)"
        print ""
        print "-h or -help                              | displays this help message."
        print ""
        print "Note on symmetry mode:"
        print "Symmetry is set separately for 'a' and 'b' directions. -s 'option' = -s 'option-option'"
        print "'mirror' only changes the z-position of the atom"
        print "'cell_center' uses the middle of the cell as an inversion point"
        print "'atom_center' uses the coordinate average of metal atoms as an inversion point"
        print "'mirror-cell_center' uses the middle of the 'a-z' plane as an inversion point"
        print "'offset' assumes that there is no inversion symmetry. The upper and lower part are shifted for correct symmetry."
        print "   -when using 'offset' mode the offset gets calculated from the input structure. When using --setOffset, use another mode"
        print "Short notations: 'mir' = 'mirror, 'cell' = 'cell_center', 'atom' = 'atom_center' and 'off' = 'offset'"


#DEFINES
def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    poscar_in = Poscar()
    poscar_out = Poscar()
    poscar_in.read(arguments.posfile)
    poscar_out.read(arguments.posfile)
    poscar_out.mirror(arguments.mode, arguments.symmetry, arguments.zplane, arguments.zrange, arguments.setOffset)
    poscar_out.move2unitcell()
    poscar_out.relabel(poscar_in,"metal")
    poscar_out.write(rewrite=True)

#EXECUTION
main(sys.argv)
