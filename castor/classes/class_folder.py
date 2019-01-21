#! /usr/bin/env python3

# 
# class_folder.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import os
import shutil

#MY PATHS

#MY CLASSES

class Folder:
    """Defines the writing to filesystem"""
    def __init__(self, root_in):
        self.root = root_in + "/"
    def make_root(self):
        os.system('mkdir ' + self.root)
    def add_file(self, source, name):
        shutil.copy2(source, self.root+name)
    def add_qsub(self, source, name, submitname):
        with open(self.root+name, 'w') as new_qsub:
            with open(source, 'r') as old_qsub:
                for line in old_qsub:
                    new_qsub.write(line.replace('#PBS -N NAME', "#PBS -N " + submitname))
    def add_poscar(self, structure_in, formatting, name):
        with open(self.root+name, 'w') as new_poscar:
            new_poscar.write(structure_in.name+"\n")
            new_poscar.write(formatting.floats([structure_in.lc]))
            new_poscar.write(formatting.floats([structure_in.vec_a.x, structure_in.vec_a.y, structure_in.vec_a.z]))
            new_poscar.write(formatting.floats([structure_in.vec_b.x, structure_in.vec_b.y, structure_in.vec_b.z]))
            new_poscar.write(formatting.floats([structure_in.vec_c.x, structure_in.vec_c.y, structure_in.vec_c.z]))
            new_poscar.write(formatting.ints(structure_in.elnr))
            new_poscar.write("Direct\n")
            for atom in structure_in.atoms:
                new_poscar.write(formatting.floats([atom.vec.x, atom.vec.y, atom.vec.z]))