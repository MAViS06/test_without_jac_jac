#!/usr/bin/python3.4

import numpy as np
import math
import scipy.linalg as sl
import random

class robot_command:

	#############################################""

	def __init__(self):
		self.mat_float_attached_points = []
		self.vec_str_names_attached_points = []
		self.vec_float_S = []
		self.vec_float_G = []
		self.vec_float_lengths = []

		self.mat_float_weights = []
		self.sc_float_lambda = .1
		self.sc_float_alpha = .6

		self.vec_float_s = []
		self.vec_float_s_star = []

		self.vec_float_s_n = []
		self.vec_float_s_nplus1 = []
		self.counter = 0 ;

	def config(self, filename):

		config_file = open(filename, 'r')

		str_mat = str()
		line = config_file.readline().rstrip('\n\r')
		nb_points = eval(line)

		for i in range(nb_points):
			line = config_file.readline().rstrip('\n\r').split(',')
			self.vec_str_names_attached_points.append(line[0])			
			str_mat += line[1] + ',' + line[2] + ',' + line[3]
			if (i < nb_points-1):
				str_mat += ';'

		self.mat_float_attached_points = np.mat(str_mat)
		config_file.close()

#		self.vec_float_lengths = np.mat('[0.0; 0.0; 0.0; 0.0]')
#		self.mat_float_weights = np.mat('[0.0, 0.0, 0.0; 0.0, 0.0, 0.0; 0.0, 0.0, 0.0; 0.0, 0.0, 0.0')
		self.vec_float_lengths = np.mat('[0.0; 0.0; 0.0; 0.0]')
		self.mat_float_weights = np.mat('[0.0, 0.0, 0.0; 0.0, 0.0, 0.0; 0.0, 0.0, 0.0 ; 0.0, 0.0, 0.0]')


	#############################################""


	def setS(self, x, y, z):
		str_coordinates = '{0}; {1}; {2}'.format(x, y, z)
		self.vec_float_S = np.mat(str_coordinates)
		self.vec_float_oldS = self.vec_float_S


	def setG(self, x, y, z):
		str_coordinates = '{0}; {1}; {2}'.format(x, y, z)
		self.vec_float_G = np.mat(str_coordinates)


	#############################################""


	def initError(self):
		self.vec_float_s = self.vec_float_S - self.vec_float_G
		self.vec_float_s_star = self.vec_float_G - self.vec_float_G
		self.vec_float_s_n = self.vec_float_G - self.vec_float_G
		self.vec_float_s_nplus1 = self.vec_float_G - self.vec_float_G


	def crossProduct(self, a, b):
		x = a[1,0] * b[2,0] - a[2,0] * b[1,0]
		y = a[2,0] * b[0,0] - a[0,0] * b[2,0]
		z = a[0,0] * b[1,0] - a[1,0] * b[0,0]
		return [x, y, z]

	def initWeights(self):
		self.mat_float_weights = np.mat('[1.0 1.0 -1.0 ; -1.0 1.0 -1.0 ; -1.0 -10 -1.0 ; 1.0 -1.0 -1.0]')* math.sqrt(3)/3
	#############################################""


	def computeLengths(self):
		for i in range(4):
			ABi = 0
			for j in range(3):
				ABi += (self.mat_float_attached_points[i,j]- self.vec_float_S[j])**2
			self.vec_float_lengths[i] = math.sqrt(ABi)


	def updateLengths(self):
		if (self.counter < 3):
			l = self.sc_float_lambda/math.exp(3-self.counter)
		else:
			l = self.sc_float_lambda

		self.vec_float_lengths += -l * self.mat_float_weights * self.vec_float_s

	#############################################""


	def computePosition(self):
		OA12 = 0
		OA22 = 0
		OA32 = 0
		for i in range(3):
			OA12 += self.mat_float_attached_points[0,i]**2
			OA22 += self.mat_float_attached_points[1,i]**2
			OA32 += self.mat_float_attached_points[2,i]**2

		Bx = ((self.vec_float_lengths[1]**2 - self.vec_float_lengths[0]**2) - (OA22 - OA12)) /(2*(self.mat_float_attached_points[0,0]-self.mat_float_attached_points[1,0])) 
		By = ((self.vec_float_lengths[2]**2 - self.vec_float_lengths[0]**2) - (OA32 - OA12)) /(2*(self.mat_float_attached_points[0,1]-self.mat_float_attached_points[2,1]))
		
		a = 1
		b = -2 * self.mat_float_attached_points[0,2]
		c = Bx**2 + By**2 + OA12 -2*(Bx*self.mat_float_attached_points[0,0] + By*self.mat_float_attached_points[0,1]) - self.vec_float_lengths[0]**2

		delta = b**2 - 4*a*c
		Bz1 = (- b - math.sqrt(delta))/(2*a)
		Bz2 = (- b + math.sqrt(delta))/(2*a)

		if Bz1 < self.mat_float_attached_points[0,2]:
			Bz = Bz1
		else:
			Bz = Bz2

		self.vec_float_S[0] = Bx ;
		self.vec_float_S[1] = By ;
		self.vec_float_S[2] = Bz ;

		self.vec_float_s = self.vec_float_S - self.vec_float_G


	#############################################""

	def updateWeights(self):

		if (self.counter > 0):
			self.vec_float_s_n = self.vec_float_s_nplus1
			self.vec_float_s_nplus1 = self.vec_float_s

			ds = self.vec_float_s_nplus1 - self.vec_float_s_n

			if (self.counter < 3):
				l = self.sc_float_lambda/math.exp(3-self.counter)
			else:
				l = self.sc_float_lambda


			dq = -l* self.mat_float_weights * self.vec_float_s			
			Ls = self.mat_float_weights

			alpha = self.sc_float_alpha

			coef_correction = alpha/(ds.transpose()*ds)
			value_correction = (dq-Ls*ds) * ds.transpose()

			mat_corr = coef_correction[0,0]*value_correction

			Ls = (1-alpha)*Ls + alpha*mat_corr

			Ls[0,:] /= np.linalg.norm(Ls[0,:])
			Ls[1,:] /= np.linalg.norm(Ls[1,:])
			Ls[2,:] /= np.linalg.norm(Ls[2,:])
			Ls[3,:] /= np.linalg.norm(Ls[2,:])

			self.mat_float_weights = Ls

			print('K après :\n', self.mat_float_weights)

		self.counter += 1

#############################################


trajectory = str()
total_error = str()
lengths = str()

rc = robot_command()
rc.config("./config2.script")

rc.setS(0.2, 1.7, 2.0)
print('Initial S :\n', rc.vec_float_S)
trajectory += '{0} {1} {2}\n'.format(rc.vec_float_S[0,0], rc.vec_float_S[1,0], rc.vec_float_S[2,0])

rc.computeLengths()
rc.computePosition()

rc.setG(0.2, 0.3, 2.0)

rc.initError()
rc.initWeights()

rc.updateLengths()
rc.computePosition()

trajectory += '{0} {1} {2}\n'.format(rc.vec_float_S[0,0], rc.vec_float_S[1,0], rc.vec_float_S[2,0])

for i in range(1024):

	rc.updateLengths()
	lengths += '{0} {1} {2} {3}\n'.format(i, rc.vec_float_lengths[0,0], rc.vec_float_lengths[1,0], rc.vec_float_lengths[2,0])

	rc.computePosition()
	trajectory += '{0} {1} {2}\n'.format(rc.vec_float_S[0,0], rc.vec_float_S[1,0], rc.vec_float_S[2,0])

	SG = rc.vec_float_G - rc.vec_float_S ;
	N = np.linalg.norm(SG)

	rc.updateWeights()

	total_error += '{0} {1}\n'.format(i, N)
	print(i, 'iterations')
	if (N < .005):
		print("Ok")
		break

print('G :\n', rc.vec_float_G)
print('Final S :\n', rc.vec_float_S) 

fic_traj = open('./trajectory.log', 'w')
fic_traj.write(trajectory)
fic_traj.close()

fic_total_error = open('./total_error.log', 'w')
fic_total_error.write(total_error)
fic_total_error.close()

fic_lengths = open('./lengths.log', 'w')
fic_lengths.write(lengths)
fic_lengths.close()