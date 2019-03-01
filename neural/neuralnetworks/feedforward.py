import time
import random
import numpy as np
random.seed(int(time.time())) # some random seed lol

"""
A feedforward neural network is an artificial neural network wherein
connections between the nodes do not form a cycle. As such, it is
different from recurrent neural networks. The feedforward neural
 network was the first and simplest type of artificial neural network devised.
"""

# sigmoid activation function
def tansig(x):
    return tanh(x)

# numpy matrix of input examples
xor_in = np.matrix([[0.0, 0.0],
         [0.0, 1.0],
         [1.0, 0.0],
         [1.0, 1.0]    ])
