import numpy as np
from matplotlib import pyplot as plt
from pauli_matrices import pauli
from pauli_matrices import pauli_projectors as projectors
default_rho = np.asarray([[1,0], [0,0]])

# The purpose of this script is to compute the PDM for a given density matrix. For now, it is programmed to use a decoherence channel. Epsilon can be passed as an argument to that function.

# The code is messy, I aim to clean it up soon.


def apply_remove_corr(rho):
	new_rho = np.zeros(rho.shape)
	np.fill_diagonal(new_rho, rho.diagonal())
	return new_rho

def apply_decoherence(rho, epsilon = 0.2):
	rho_diag = apply_remove_corr(rho)
	new_rho = (1 - epsilon) * rho
	new_rho = new_rho + epsilon * rho_diag
	return new_rho

def compute_expectations(rho, i, k, epsilon):
	rho_original = rho;
	expectation = 0
	operator = pauli[i]
	for j in range(2):
		rho = rho_original
		projector = projectors[i,j]
		probability = np.trace(np.dot(rho, projector))
		if probability.real == 0:
			continue
		rho = np.dot(rho, projector)
		rho = np.dot(projector, rho)
		rho = rho/float(probability)
		eigenvalue = (-1)**j
		if i == 3:
			eigenvalue = 1
		rho = apply_decoherence(rho, epsilon)
		operator2 = pauli[k]
		rho_original_second = rho
		for j2 in range(2):
			rho = rho_original_second
			projector2 = projectors[k,j2]
			probability2 = np.trace(np.dot(rho,projector2))
			if probability2 == 0:
				continue
			else:
				rho = np.dot(rho, projector2)
				rho = np.dot(projector2, rho)
				rho = rho/float(probability2)
			eigenvalue2 = (-1)**j2
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
			#print expectation
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


