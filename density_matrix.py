import numpy as np
from pauli_matrices import pauli
from pauli_matrices import pauli_projectors as projectors
default_rho = np.asarray([[1,0], [0,0]])

# The purpose of this script is to compute the PDM for a given density matrix. For now, it is programmed to use a decoherence channel. Epsilon can be passed as an argument to that function.

# The code is messy, I aim to clean it up soon.


def apply_remove_corr(rho):
	new_rho = np.zeros(rho.shape)
	np.fill_diagonal(new_rho, rho.diagonal())
	return new_rho

def apply_decoherence(rho, epsilon = 0):
	rho_diag = apply_remove_corr(rho)
	new_rho = (1 - epsilon) * rho
	new_rho = new_rho + epsilon * rho_diag
	return new_rho

def compute_expectations(rho, i, k):
	rho_original = rho;
	counter = 1
	expectation = 0
	if i != 3 and k != 3:
		#print 'Here we start'
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
			rho = apply_decoherence(rho)
			#rho = apply_remove_corr(rho)
			#Add unitary evolution here
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
				print j2
				counter = counter + 1
	elif i == 3 and k != 3:
		eigenvalue = 1
		probability = 1
		rho = apply_decoherence(rho)
		#rho = apply_remove_corr(rho)
		#Second measurement starts here
		operator2 = pauli[k]
		for j2 in range(2):
			rho = rho_original
			projector2 = projectors[k,j2]
			probability2 = np.trace(np.dot(rho,projector2))
			if probability2 == 0:
				continue
			rho = np.dot(rho, projector2)
			rho = np.dot(projector2, rho)
			rho = rho/float(probability2)
			eigenvalue2 = (-1)**j2
			expectation = expectation + probability*probability2*eigenvalue*eigenvalue2
	elif k == 3 and i != 3:
		operator = pauli[i]
		for j in range(2):
			rho = rho_original
			projector = projectors[i,j]
			probability = np.trace(rho * projector)
			if probability == 0:
				continue
			rho = np.dot(rho, projector)
			rho = np.dot(projector, rho)
			rho = rho/float(probability)
			eigenvalue = (-1)**j
			#rho = apply_remove_corr(rho)
			rho = apply_decoherence(rho)
			eigenvalue2 = 1
			probability2 = 1
			expectation = expectation + probability*probability2*eigenvalue*eigenvalue2
	else:
		expectation = 1
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

def compute_pdm(rho, channel_unitary = np.asarray([[1,0],[0,1]])):
	pdm = np.zeros([4,4])
	for i in range(4):
		for jj in range(4):
			expectation = compute_expectations(rho, i, jj)
			#print expectation
			basis_matrix = np.kron(pauli[i],pauli[jj])
			temp = expectation * basis_matrix
			pdm = pdm + temp
	#Normalising the pdm
	pdm = pdm/np.trace(pdm)
	return pdm
