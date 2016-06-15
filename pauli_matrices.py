import numpy as np

# This script just contains useful matrices to be used later.

X = np.asarray([[0,1], [1,0]])
Y = np.asarray([[0, -1j], [1j, 0]])
Z = np.asarray([[1,0], [0,-1]])
I = np.asarray([[1,0], [0,1]])

pauli = [X,Y,Z,I]

xp1 = 0.5*np.asarray([[1,1], [1,1]])
xp2 = 0.5*np.asarray([[1,-1], [-1,1]])

yp1 = 0.5*np.asarray([[1, -1j], [1j, 1]])
yp2 = 0.5*np.asarray([[1, 1j], [-1j, 1]])

zp1 = np.asarray([[1, 0], [0, 0]])
zp2 = np.asarray([[0, 0], [0, 1]])

ip1 = np.identity(2)
ip2 = np.zeros([2,2])

pauli_projectors = np.zeros([4,2,2,2])
pauli_projectors = pauli_projectors + 0j
pauli_projectors[0,0] = xp1
pauli_projectors[0,1] = xp2
pauli_projectors[1,0] = yp1
pauli_projectors[1,1] = yp2
pauli_projectors[2,0] = zp1
pauli_projectors[2,1] = zp2
pauli_projectors[3,0] = ip1
pauli_projectors[3,1] = ip2

