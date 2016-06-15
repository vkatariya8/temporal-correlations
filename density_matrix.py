import numpy as np
import math
from matplotlib import pyplot as plt
from pauli_matrices import pauli
from pauli_matrices import pauli_projectors as projectors
default_rho = np.asarray([[1,0], [0,0]])

# The purpose of this script is to compute the PDM for a given density matrix. For now, it is programmed to use a decoherence channel. Epsilon can be passed as an argument to that function.

# The code is messy, I aim to clean it up soon.


def operator_sum_apply(rho, operators):
	newrho = np.zeros(rho.shape)
	for i in range(len(operators)):
		temprho = np.dot(operators[i], np.dot(rho, np.asmatrix(operators[i]).H))
		newrho = newrho + temprho
	return newrho


def apply_channel(rho, option = 1, param1 = 0, param2 = 0):
	'''
	Options 1: decoherence channel, 2: dephasing, 3: depolarising
	1: param1 = epsilon, param2 unused
	'''
	if option == 1:
		rho = apply_decoherence(rho, param1)
	elif option == 2:
		rho = apply_dephasing(rho, param1)
	elif option == 3:
		rho = apply_depolarising(rho, param1)
	return rho

def apply_depolarising(rho, p):
	rho = p * np.identity(2) / 2 + (1 - p) * rho
	return rho


def apply_dephasing(rho, param):
	operators = np.zeros([2,2,2])
	operators[0] = np.asarray([[1,0],[0,math.sqrt(1 - param)]])
	operators[1] = np.asarray([[0,0],[0,math.sqrt(param)]])
	rho = operator_sum_apply(rho, operators)
	return rho

def apply_remove_corr(rho):
	new_rho = np.zeros(rho.shape)
	np.fill_diagonal(new_rho, rho.diagonal())
	return new_rho

def apply_decoherence(rho, epsilon = 0.2):
	rho_diag = apply_remove_corr(rho)
	new_rho = (1 - epsilon) * rho
	new_rho = new_rho + epsilon * rho_diag
	return new_rho

def perform_measurement(rho, i, j):
	projector = projectors[i,j]
	probability = np.trace(np.dot(rho, projector))
	if probability.real == 0:
		return (0,1,rho)
	rho = np.dot(projector, np.dot(rho, projector))
	rho = rho/float(probability)
	eigenvalue = (-1)**j
	if i == 3:
		eigenvalue = 1
	return (probability, eigenvalue, rho)

def compute_expectations(rho, i, k, epsilon):
	rho_original = rho;
	expectation = 0
	for j in range(2):
		rho = rho_original
		(probability, eigenvalue, rho) = perform_measurement(rho, i, j)
		rho = apply_channel(rho, 2, epsilon)
		rho_original_second = rho
		for j2 in range(2):
			rho = rho_original_second
			(probability2, eigenvalue2, rho) = perform_measurement(rho, k, j2)
			expectation = expectation + probability.real*probability2.real*eigenvalue*eigenvalue2
	return expectation

def define_density_operator(a = 1,b = 0,c = 0,d = 0):
	if (a + d) != 1:
		print 'Trace not equal to 1, beware of false results'
		rho = default_rho
	if b != np.conjugate(c):
		print 'Not hermitian, beware of false results'
		rho = default_rho
	rho = np.asarray([[a,b], [c,d]])
	return rho

def compute_pdm(rho, epsilon):
	pdm = np.zeros([4,4])
	for i in range(4):
		for jj in range(4):
			expectation = compute_expectations(rho, i, jj, epsilon)
			basis_matrix = np.kron(pauli[i],pauli[jj])
			temp = expectation * basis_matrix
			pdm = pdm + temp
	#Normalising the pdm
	pdm = pdm/np.trace(pdm)
	return pdm

def pdm_analysis(pdm):
	eigen_values = np.linalg.eigvals(pdm)
	smallest = min(eigen_values)
	abs_values = abs(eigen_values)
	smallest_mod = min(abs_values)
	return (smallest, smallest_mod)

def iterate_over_epsilon(rho):
	lower_epsilon = 0
	upper_epsilon = 1
	step_size = 0.01
	iterations = (upper_epsilon - lower_epsilon)/step_size
	iterations = int(iterations)
	small_list = list()
	abs_small_list = list()
	eigvalues = np.zeros([4, iterations])
	for i in range(iterations):
		epsilon = lower_epsilon + i*upper_epsilon*step_size
		pdm = compute_pdm(rho, epsilon)
		#(small, smallmod) = pdm_analysis(pdm)
		#print small
		#small_list.append(small)
		eigvals = np.linalg.eigvals(pdm)
		for j in range(4):
			temp = eigvals[j]
			eigvalues[j][i] = temp
	return eigvalues

def plot_eigvalues(eigvalues):
	for i in range(len(eigvalues)):
		plt.plot(eigvalues[i])
	plt.show()


