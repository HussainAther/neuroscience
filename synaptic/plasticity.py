import numpy as np

from sklearn import make_blobs

"""
Synaptic plasticity with Hebbian Learning. Using neurons of the brain as processing units
using cable theory as extensions, we can simulate models of neurons firing using action potentials
that deliber chemical siganls (neurotransmitters) as ion channels open and close to allow
ions to flow.

"Neurons that fire together wire together."
With Hebbian theory, we can calculate the synaptic weight from a neuron if both pre-synaptic
and post-synaptic units behave in the same way (whether they"re firing or in steady state).
Hebbian learning uses principal component analysis by extracting the first principal component.
Hebb"s rule can be expressed as:

Δw_i = α x_i y

in which we use w as the weight, α as the learning rate, x as each value, and y as the linear combination of w and
x for a particular input value. This rule is unstable because when expressed using vectors (Δw = αC . w), when
we solve for the eigenvalues of C (a positive semidefinite matrix), we get eigenvectors that increase
exponentially with t. This will cause an overflow with x and y values greater than 1.

Finnish computer scientist Erkki Oja proposed an alternative. Like Hebb"s rule, we use the vector w that converges
to the C eigenvector, but, in this case, to a finite (small) number.

Δw_i = α y(x̄_i - w̅_iy)
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

# Oja"s rule
woja = np.random.normal(scale=.25, size=(2,1)) # get the weights for a normal distribution
prevwoja = np.ones((2,1)) # previous Oja weights to keep track as we go along

a = .0001 # alpha: learning rate
eps = 1e-8 # epsilon: tolerance

while np.linalg.norm(prevwoja - woja) > eps:
    # perform the iteration to calculate the Ys values of the Oja weights. 
    prevwoja = woja.copy()
    Ys = np.dot(Xs, woja) # get the dot-product with the Oja weights
    woja += a * np.sum(Ys*Xs - np.square(Ys)*woja.T, axis=0).reshape((2,1)) # perform the calculation and add to overall Oja weights

"""
We may extend Oja"s rule to multi-output networks by the Sanger (Sanger"s) rule
which is also known as the Generalized Hebbian Algorithm.

Δw_ij = α (y_i x_j - y_i summation from k=1 to i of w_kj y_k)

such that the update rule is now

Δw =  α (y̅ . x̄^T - tril(y̅ x y̅^T) . w̅)

in which tril is a function that returns the lower triangle of a square matrix.
"""

n = 3000 # number of iterations
t = 0 # initial time

wsanger = np.random.normal(scale=.01, size(2,2)) # Sanger"s rule
prevw = np.ones((2,2)) # keep track of previous values

for i in range(n):
    prevw = wsanger.copy() # store the previous value for the next iteration. 
    dw = np.zeros((2,2)) # delta w
    t += 1
    for j in range(Xs.shape[0]):
        Ysj = np.dot(sanger, Xs[j]).reshape((2,1))
        tril = np.tril(np.dot(Ysj, Ysj.T)) # lower triangle of square matrix
        dw += np.dot(Ysj, Xs[j].reshape(1,2))) - np.dot(tril, wsanger)

    wsanger += (alpha/t) * dw
    wsanger /= np.linalg.norm(wsanger, axis=1).reshape((2,1))
