import numpy as np
import matplotlib.pyplot as plt
from density_matrix import compute_expectations
from pauli_matrices import pauli


def compute_pdm(rho1 = 0.5*np.identity(2), rho2 = 0.5*np.identity(2)):
	pdm = np.zeros([16,16])
	pdm = pdm + 0j
	for i in range(4):
		for j in range(4):
			for k in range(4):
				for l in range(4):
					expectation1 = compute_expectations(rho1, i, j, 0)
					expectation2 = compute_expectations(rho2, k, l, 0)
					expectation = expectation1*expectation2
					basis_1 = np.kron(pauli[i], pauli[j])
					basis_2 = np.kron(pauli[k], pauli[l])
					basis_m = np.kron(basis_1, basis_2)
					temp = expectation*basis_m
					pdm += temp
	pdm /= np.trace(pdm)
	return pdm
