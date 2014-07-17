import numpy as np
import math
import scipy.linalg as sl

class robot_command:

	#############################################""

	def __init__(self):
		self.mat_float_attached_points = []
		self.vec_str_names_attached_points = []
		self.vec_float_S = []
		self.vec_float_G = []
		self.vec_float_lengths = []

		self.mat_float_weights = []
		self.vec_float_weights = []
		self.sc_float_delta_lambda = .1

		self.vec_float_s1 = []
		self.vec_float_s2 = []
		self.vec_float_s_star = []

		self.sc_float_alpha_weight_correction = .01

		self.vec_float_correction = []

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

		self.vec_float_lengths = np.mat('[0.0; 0.0; 0.0]')
		self.mat_float_weights = np.mat('[0.0, 0.0, 0.0; 0.0, 0.0, 0.0; 0.0, 0.0, 0.0]')
		self.vec_float_weights = np.mat('[0.0; 0.0; 0.0]')
		self.vec_float_correction = np.mat('[0.0; 0.0; 0.0]')


	#############################################""


	def setS(self, x, y, z):
		str_coordinates = '{0}; {1}; {2}'.format(x, y, z)
		self.vec_float_S = np.mat(str_coordinates)


	def setG(self, x, y, z):
		str_coordinates = '{0}; {1}; {2}'.format(x, y, z)
		self.vec_float_G = np.mat(str_coordinates)


	#############################################""


	def initError(self):
		self.vec_float_s1 = self.vec_float_G - self.vec_float_S
		self.vec_float_s2 = self.vec_float_G - self.vec_float_S
		self.vec_float_s_star = self.vec_float_G - self.vec_float_G

	def initWeights(self):
		for i in range(3):
			for j in range(3):
				self.mat_float_weights[i,j] = np.copysign(1.0/3.0, (self.vec_float_G[j] - self.vec_float_S[j]) * (self.vec_float_G[j] - self.mat_float_attached_points[i,j]))
				self.vec_float_weights[i] += self.mat_float_weights[i,j]


	#############################################""


	def computeLengths(self):
		for i in range(3):
			ABi = 0
			for j in range(3):
				ABi += (self.mat_float_attached_points[i,j]- self.vec_float_S[j])**2
			self.vec_float_lengths[i] = math.sqrt(ABi)


	def updateLengths(self):
		self.vec_float_lengths += self.sc_float_delta_lambda * self.vec_float_weights


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

		self.vec_float_s1 = self.vec_float_s2
		self.vec_float_s2 = self.vec_float_G - self.vec_float_S


	#############################################""


	def updateWeights(self):
		for i in range(3):
			self.vec_float_correction[i] = (self.vec_float_s2[i] - self.vec_float_s_star[i])/(self.vec_float_s1[i] - self.vec_float_s_star[i])

		for i in range(3):
				self.mat_float_weights[:,i] = self.mat_float_weights[:,i].dot(self.vec_float_correction[i])

		for i in range(3):
			self.vec_float_weights[i] = 0
			for j in range(3):
				self.vec_float_weights[i] += self.mat_float_weights[i,j]

#############################################


trajectory = str()
total_error = str()
lengths = str()

rc = robot_command()
rc.config("./config.script")

rc.setS(0.8, 1.0, 2.20)
trajectory += '{0} {1} {2}\n'.format(rc.vec_float_S[0,0], rc.vec_float_S[1,0], rc.vec_float_S[2,0])

rc.computeLengths()
rc.computePosition()
print('Initial S :\n', rc.vec_float_S)
rc.setG(1.1, 0.3, 1.20)

rc.initError()
rc.initWeights()

rc.updateLengths()
rc.computePosition()
trajectory += '{0} {1} {2}\n'.format(rc.vec_float_S[0,0], rc.vec_float_S[1,0], rc.vec_float_S[2,0])

for i in range(540):
	rc.updateWeights()
	rc.updateLengths()
	lengths += '{0} {1} {2} {3}\n'.format(i, rc.vec_float_lengths[0,0], rc.vec_float_lengths[1,0], rc.vec_float_lengths[2,0])

	rc.computePosition()
	trajectory += '{0} {1} {2}\n'.format(rc.vec_float_S[0,0], rc.vec_float_S[1,0], rc.vec_float_S[2,0])

	SG = rc.vec_float_G - rc.vec_float_S ;
	N = math.sqrt(SG[0]**2 + SG[1]**2 + SG[2]**2)

	total_error += '{0} {1}\n'.format(i, N)

	if (N < .005):
		print("Ok")
		print(i, 'iterations')
		print('Final S :\n', rc.vec_float_S)
		break

print('G :\n', rc.vec_float_G)


fic_traj = open('./trajectoire.log', 'w')
fic_traj.write(trajectory)
fic_traj.close()

fic_total_error = open('./total_error.log', 'w')
fic_total_error.write(total_error)
fic_total_error.close()

fic_lengths = open('./lengths.log', 'w')
fic_lengths.write(lengths)
fic_lengths.close()