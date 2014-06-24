#! /usr/bin/env python

# 
# class_formatting.py
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

class Formatting:
    """Defines settings for Potcar formatting compatible with vasp"""
    def __init__(self, decimals_in, spaces_in):
        self.decimals = decimals_in
        self.spaces = spaces_in
    def floats(self, floatlist):
        string = ""
        for number in floatlist:
            temp = '{:.{width}f}'.format(number, width=self.decimals)
            string += '{:>{width}}'.format(temp, width=self.decimals+self.spaces+2)
        string += "\n"
        return string
    def ints(self, intlist):
        firstnum = 1
        string = ""
        for number in intlist:
            temp = '{:.{width}f}'.format(number, width=self.decimals)
            string += '{:>{width}}'.format(temp, width=self.decimals+self.spaces+2)
            if(firstnum==1):
                string = str(number)
                firstnum = 0
            else:
                for i in range(0, self.spaces):
                    string += " "
                string += str(number)
        string += "\n"
        return string