import numpy as np

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
"""


