import numpy as np

# Construct W, the network weight matrix.
W = np.ones((5,5))
W = W / 10.
np.fill_diagonal(W, 0.6)

# Construct u, the static input vector.
u = np.zeros(5)
u[0] = 0.6
u[1] = 0.5
u[2] = 0.6
u[3] = 0.2
u[4] = 0.1

# Construct M, the recurrent weight matrix.
M = np.zeros((5,5))
np.fill_diagonal(M, -0.25)
for i in range(3):
    M[2+i][i] = 0.25
    M[i][2+i] = 0.25
for i in range(2):
    M[3+i][i] = 0.25
    M[i][3+i] = 0.25

# We need to matrix multiply W and u together to get h
# NOTE: cannot use W * u, that"s going to do a scalar multiply
# it"s element wise otherwise
h = W.dot(u)

# Ok then the big deal is:
#                               h dot e_i
# v_ss = sum_(over all eigens) ------------ e_i
#                               1 - lambda_i

eigs = np.linalg.eig(M)

eigenvalues = eigs[0]
eigenvectors = eigs[1]

v_ss = np.zeros(5)
for i in range(5):
    v_ss += (np.dot(h,eigenvectors[:, i]))/((1.0-eigenvalues[i])) * eigenvectors[:,i]
