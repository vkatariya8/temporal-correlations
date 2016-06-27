import numpy as np
from tabulate import tabulate
import math
from matplotlib import pyplot as plt
from pauli_matrices import pauli
from density_matrix import apply_decoherence, perform_measurement, plot_eigvalues

correlations = [1,-1,1]

def compute_rho_from_eigenstate(i, ii):
	if i == 0:
		if ii == 0:
			rho = 0.5*np.asarray([[1,1], [1,1]])
		elif ii == 1:
			rho = 0.5*np.asarray([[1,-1],[-1,1]])
	elif i == 1:
		if ii == 0:
			rho = 0.5*np.asarray([[1,-1j],[1j,1]])
		elif ii == 1:
			rho = 0.5*np.asarray([[1,1j], [-1j, 1]])
	elif i == 2:
		if ii == 0:
			rho = np.asarray([[1,0], [0,0]])
		elif ii == 1:
			rho = np.asarray([[0,0], [0,1]])
	elif i == 3:
		rho = 0.5*np.identity(2)
	return rho

def triple_expectation(i, j, k, e):
	expectation = 0
	for ii in range(2):
		eigenvalue = (-1)**ii
		if i == 3:
			eigenvalue = 1
		probability = 0.5
		for jj in range(2):
			eigenvalue2 = (-1)**jj
			probability2 = 0.5
			if j == 3:
				eigenvalue2 = 1
			elif i == j:
				eigenvalue2 = correlations[i]*eigenvalue
			# what state is the first qubit in now? has anything changed from the fist measurement? no, right
			# have to now apply the channel - have to reuse the code written for density matrix
			rho = compute_rho_from_eigenstate(i,ii)
			rho = apply_decoherence(rho, e)
			# for kk in range(2):
			# 	eigenvalue3 = (-1)**kk
			# 	if k == 3:
			# 		eigenvalue3 = 1
			# 	probability3 = 0.5
			# 	if i == k:
			# 		eigenvalue3 = eigenvalue
			# 	tempe = probability*probability2*probability3*eigenvalue*eigenvalue2*eigenvalue3
			# 	expectation += tempe
			rho_original = rho
			for kk in range(2):
				rho = rho_original
				(probability3, eigenvalue3, rho) = perform_measurement(rho, k, kk)
				tempe = probability*probability2*probability3*eigenvalue*eigenvalue2*eigenvalue3
				#tempe = probability*probability3*eigenvalue*eigenvalue3
				#print probability, probability3, eigenvalue, eigenvalue3
				expectation += tempe
	return expectation

def compute_pdm(epsilon):
	pdm = np.zeros([8,8])
	for i in range(4):
		for j in range(4):
			for k in range(4):
				exp = triple_expectation(i,j,k, epsilon)
				basis_m = np.kron(pauli[i], pauli[j])
				basis_m = np.kron(basis_m, pauli[k])
				#basis_m = np.kron(pauli[i], pauli[k])
				pdm = pdm +  exp*basis_m
	#pdm = pdm/4.0
	#pdm = pdm/4.0
	pdm = pdm/np.trace(pdm)
	return pdm

def iterate_epsilon():
	lower_epsilon = 0
	upper_epsilon = 1
	step_size = 0.001
	iterations = int((upper_epsilon - lower_epsilon)/step_size)
	epsilons = np.linspace(lower_epsilon, upper_epsilon, iterations)
	eigenvalues = np.zeros([8, iterations])
	trace_norms = np.zeros(iterations)
	for i in range(iterations):
		epsilon = lower_epsilon + step_size * i
		p = compute_pdm(epsilon)
		eigvals = np.linalg.eigvalsh(p)
		trace_norm = sum(np.abs(eigvals))
		trace_norms[i] = trace_norm
		for j in range(8):
			eigenvalues[j,i] = eigvals[j]
	return (eigenvalues, trace_norms, epsilons)

# def iterate_coefficients():
# 	la = 0
# 	ua = 0.4
# 	global coefficients
# 	step_size = 0.01
# 	iterations = (ua - la)/step_size
# 	iterations = int(iterations)
# 	for i in range(iterations):
# 		alpha = la + i*step_size
# 		beta = math.sqrt(1 - alpha**2)
# 		coefficients = [alpha, beta]
# 		(p,t) = iterate_epsilon()
# 		plt.plot(t)
# 		#plot_eigvalues(p)
# 	plt.show()

