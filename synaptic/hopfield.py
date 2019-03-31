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

We can define hte energy of a Hopfield network as :

E = (-1/2)*(summation of i)(summation of j) of w_ij*x_i*x_j + (summation of i) of θ_i * x_i

We're going to use Hopfield networks to improve the pattern performance of a few 5x5 grids patterns. 
"""

p = 4
