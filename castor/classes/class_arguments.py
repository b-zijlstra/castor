#! /usr/bin/env python3

# 
# class_arguments.py
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

class Arguments:
    """Defines a set of arguments"""
    def __init__(self, defaults_in, options_in, help_in):
        self.defaults = defaults_in
        self.options = options_in
        self.help = help_in
    def setup(self, arg_in):
        self.readmode = None
    def help(self):
        pass