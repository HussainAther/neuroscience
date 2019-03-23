import numpy as np

from sklearn import make_blobs
"""
Synaptic plasticity with Hebbian Learning. Using neurons of the brain as processing units
using cable theory as extensions, we can simulate models of neurons firing using action potentials
that deliber chemical siganls (neurotransmitters) as ion channels open and close to allow
ions to flow.

"Neurons that fire together wire together."
With Hebbian theory, we can calculate the synaptic weight from a neuron if both pre-synaptic
and post-synaptic units behave in the same way (whether they're firing or in steady state).
Hebbian learning uses principal component analysis by extracting the first principal component.
Hebb's rule can be expressed as:

Δw_i = α x_i y


in which we use w as the weight, α as the learning rate, x as each value, and y as the linear combination of w and
x for a particular input value. This rule is unstable because when expressed using vectors (Δw = αC . w), when
we solve for the eigenvalues of C (a positive semidefinite matrix), we get eigenvectors that increase
exponentially with t. This will cause an overflow with x and y values greater than 1.

Finnish computer scientist Erkki Oja proposed an alternative. Like Hebb's rule, we use the vector w that converges
to the C eigenvector, but, in this case, to a finite (small) number.

Δw_i = α y(x_i - w_iy)
"""

# create the dataset
X, y =  make_blobs(n_samples=500, centers=2, cluster_std=5.0, random_state=1000)

def scale(x):
    """
    Transform dataset x such that it will have a mean value 0 and
    standard deviation 1.
    """
    a = np.mean(x) # find the mean/average of the dataset
    s = np.std(x) # find the standard deviation of the dataset
    r = [] # output result list
    for i in x:
        r.append((i-a)/s)
    return r

# scale the dataset
Xs = scale(X)

# compute the covariance, eigenvalues, and eigenvectors
q = np.cov(Xs.T)
u, v = np.linalg.eig(q)

# Oja's rule

