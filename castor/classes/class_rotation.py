#! /usr/bin/env python3

# 
# class_rotation.py
# 
# Author: Bart Zijlstra
# 
# (C) Copyright 2014 Inorganic Materials Chemistry
# 
# 

#IMPORTS
import sys
import math

#MY PATHS

#MY CLASSES

class Matrix:
	def __init__(self, val):
		self.v = []
		for i in range(0,9):
			self.v.append(val)
	def g(self, row, col):
		return self.v[(row-1)*3 + col-1]
	def s(self, row, col, val):
		self.v[(row-1)*3 + col-1] = float(val)
	def write(self):
		for i in range(0,3):
				print(str(self.g(i+1,1))+"\t"+str(self.g(i+1,2))+"\t"+str(self.g(i+1,3)))
	def setRotationZ(self, angle):
		self.s(1,1,math.cos(angle));
		self.s(1,2,-math.sin(angle));
		self.s(1,3,0);
		self.s(2,1,math.sin(angle));
		self.s(2,2,math.cos(angle));
		self.s(2,3,0);
		self.s(3,1,0);
		self.s(3,2,0);
		self.s(3,3,1);

class AxisAngle:
	"""Takes as input an axis and an angle"""
	def __init__(self, vxx, vyy, vzz, aa):
		aa *= 2.0 * math.pi / 360.0
		n = 1.0 / math.sqrt(vxx * vxx + vyy * vyy + vzz * vzz)
		self.x = n * vxx * math.sin(aa / 2.0);
		self.y = n * vyy * math.sin(aa / 2.0);
		self.z = n * vzz * math.sin(aa / 2.0);
		self.w = math.cos(aa / 2.0)
	def toMatrix(self):
		# convert the axis-angle to quaternion and convert
		# the quaternion to a rotation matrix
		qx = self.x
		qy = self.y
		qz = self.z
		qw = self.w
		n = 1.0 / math.sqrt(qx*qx+qy*qy+qz*qz+qw*qw)
		qx *= n
		qy *= n
		qz *= n
		qw *= n
		# construct the rotation matrix
		mat = Matrix(0)
		mat.s(1,1,1.0 - 2.0*qy*qy - 2.0*qz*qz)
		mat.s(1,2,2.0*qx*qy - 2.0*qz*qw)
		mat.s(1,3,2.0*qx*qz + 2.0*qw*qy)
		mat.s(2,1,2.0*qx*qy + 2.0*qz*qw)
		mat.s(2,2,1.0 - 2.0*qx*qx - 2.0*qz*qz)
		mat.s(2,3,2.0*qy*qz - 2.0*qx*qw)
		mat.s(3,1,2.0*qx*qz - 2.0*qy*qw)
		mat.s(3,2,2.0*qy*qz + 2.0*qx*qw)
		mat.s(3,3,1.0 - 2.0*qx*qx - 2.0*qy*qy)
		return mat