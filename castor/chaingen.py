#! /usr/bin/env python3

# 
# chaingen.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2017 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import re

#MY PATHS

#MY CLASSES

#CLASSES
class Arguments:
    """Defines a set of arguments"""
    def __init__(self):
        self.mkm = "input.mkm"
        self.length = 50
    def setup(self, arg_in):
        readmode = None
        for i in range(1,len(arg_in)):
            if(readmode=="mkm"):
                self.mkm = arg_in[i]
                readmode = None
                continue
            if(readmode=="length"):
                self.length = int(arg_in[i])
                readmode = None
                continue
            if(sys.argv[i] == "-h" or sys.argv[i] == "--help"):
                self.help()
                sys.exit()
            elif(sys.argv[i] == "-i" or sys.argv[i] == "--mkm"):
                readmode = "mkm"
                continue
            elif(sys.argv[i] == "-l" or sys.argv[i] == "--length"):
                readmode = "length"
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
        print("Use: chaingen.py <options>")
        print("Options:")
        print("-i or --mkm <input name>         | Example: $chaingen.py -i base.mkm (Default = input.mkm)")
        print("-l or --length <length>          | Example: $chaingen.py -l 10 (Default = 50)")
        print("")
        print("-h or --help                     | displays this help message")

#DEFINES

def makemkm(mkm_in, length_in):
    with open(mkm_in, 'r') as mkm:
        regex_begin='^[ \t]*@BEGIN[ \t]*;[ \t]*([0-9]+)[ \t]*;[ \t]*([0-9]+)[ \t]*$'
        regex_swap='^.*@.*$'
        regex_swap_p='(^.*)@p([0-9]+)(.*$)'
        regex_swap_m='(^.*)@m([0-9]+)(.*$)'
        regex_swap_0='(^.*)@(.*$)'
        regex_end='^[ \t]*@END[ \t]*$'
        lines = mkm.readlines()
        readmode = 0
        start = 0
        end = 0
        swap = []
        for line in lines:
            line = line.rstrip("\n")
            if(readmode == 0): # Check for @BEGIN
                match = re.search(regex_begin,line)
                if(match): # start buffer
                    start = int(match.group(1))
                    end = int(match.group(2))
                    if(end==0):
                        end = length_in
                    readmode = 1
                    # print(line)
                    continue
                print(line) # just print the unmodified line
                continue
            elif(readmode == 1): # Check for substitutions or @END
                match = re.search(regex_end,line)
                if(match): # stop checking for substitutions and initiate writing
                    readmode = 0
                    for i in range(start,end+1):
                        for singleline in swap:
                            match = re.search(regex_swap,singleline)
                            if(match): # substitute line
                                subline = singleline
                                match = re.search(regex_swap_p,subline)
                                while(match):
                                    #substitute line
                                    left = match.group(1)
                                    num = int(match.group(2))
                                    right =  match.group(3)
                                    subline = left + str(i + num) + right
                                    match = re.search(regex_swap_p,subline)
                                match = re.search(regex_swap_m,subline)
                                while(match):
                                    #substitute line
                                    left = match.group(1)
                                    num = int(match.group(2))
                                    right =  match.group(3)
                                    subline = left + str(i - num) + right
                                    match = re.search(regex_swap_m,subline)
                                match = re.search(regex_swap_0,subline)
                                while(match):
                                    #substitute line
                                    left = match.group(1)
                                    right =  match.group(2)
                                    subline = left + str(i) + right
                                    match = re.search(regex_swap_0,subline)
                                print(subline)
                                continue
                            print(singleline)
                    swap = []
                    continue
                swap.append(line)
                continue
        if(readmode == 1):
            print("# WARNING: End of file reached before @END")

def main(arg_in):
    arguments = Arguments()
    arguments.setup(arg_in)
    makemkm(arguments.mkm,arguments.length)

#EXECUTION
main(sys.argv)