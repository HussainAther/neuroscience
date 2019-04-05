import numpy as np
import matplotlib.pyplot as plt

"""
John Hopfield combined a generalized form of the Hebb rule with simple activation
dynamics to construct an energy function in which the stored memory states
were local minima of the energy function. It allowed Hebbian learning to lead
to attrcator memories and strengthened the bridge between memory networks
and certain branches of statistical mechanics.

Hopfield networks are recurrent neural networks with bipolar threshold neurons.
They're used in associative memory. The synaptic weights have the following conditions:

w_ij = wji ∀ i,j ∈ {1,...,N}
w_ii = 0 ∀ i ∈ {1,...,N}

For N neurons, the activation function for each neuron is:

f_A(x) = 1 if summation from j=1 to N of w_ij*x_j > θ_i 
        -1 otherwise

We can define the energy of a Hopfield network as :

E = (-1/2)*(summation of i)(summation of j) of w_ij*x_i*x_j + (summation of i) of θ_i * x_i

We're going to use Hopfield networks to improve the pattern performance of a few 5x5 binary grids patterns.
By binary, we're using 1's and -1's in accordance with the activation functions. 
"""

p = 5 # number of patterns
w = 5 # width
l = 5 # length
iter = 20 # number of iterations

x = np.zeros(p, w, *l)) # data of binary grid pattern 

# make some patterns
x[0] = [-1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, -1, -1, -1, -1]
x[1] = [-1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, -1, 1, -1]
x[2] = [-1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1]
x[3] = [1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1]
x[4] = [1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, 1]

# plot 'em
fig, ax = plt.subplots(1, p, figsize=(15,10))

for i in range(p):
    ax[i].matshow(x[i].reshape((w, l)), cmap="gray") # show the axis w/color map cmap
    ax[i].set_xticks([])
    ax[i].set_yticks([])

plt.show()

# train
w = np.zeros((w * l, w*l))

for i in range(w*l):
    for j in range(w*l):
        if i == j or w[i, j] != 0:
            continue
        w = 0 # weight for the synapse
        for n in range(p):
            w += x[n, i] * x[n, j] # add each value to the weight as the summation dictates
        w[i, j] = w/x.shape[0] # normalize with resepct to the size of our grid 
        w[j, i] = w[i, j] # as expressed in the equation for the synaptic weights 

# try to recover the pattern
x_test = np.array([1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 1, 1, 1])

# score the test
for i in range(iter):
    for j in range(w * l):
        x_test = 1 if np.dot(w[i], x_test) > 0 else -1

# show the recovered pattern
fig, ax = plt.subplots(1, 2, figsize=(15,10))

# plot 'em
ax[1].matshow(x_test.reshape(h, l), cmap="gray")
ax[1].set_title("Recovered pattern")
ax[1].set_xticks([])
ax[1].set_yticks([])

plt.show()

"""
We can use a stochastic process on a probability distribution using the sigmoid function.
Each neuron is activated according to a probability function.
"""

def sigmoid(x):
    """
    Return the value of the sigmoid function.
    """
    return 1 / (1 + np.exp(-x))

def prob(N, i, w, x):
    """
    Each neuron has a threshold set to null that can be activated for N number of neurons at some row i with
    some weight matrix w and x values (activation state).
    """
    summ = 0
    for j in range(N):
          summ += w[i][j] * x[j] # sum each weight with each x value (as described above).
    return sigmoid(summ) 
