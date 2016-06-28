from qutip import *

rho = tensor(fock_dm(2,0), fock_dm(2, 0), fock_dm(2,0), fock_dm(2,0))
sa = scipy.sparse.csr_matrix(p, dtype = np.complex128)
rho.data = sa


# table to latex

print " \\\\\n".join([" & ".join(map(str,line)) for line in a])

print " \\\\\n".join([" & ".join(map('{0:.3f}'.format, line)) for line in a])
