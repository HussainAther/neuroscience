import time
import random
import numpy as np
import matplotlib.pyplot as plt
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

# train the matrix with these outputs
xor_out = matrix([[0.0],
          [1.0],
          [1.0],
          [0.0]    ])

# initialize weights and biases to small random values
sigw = 0.5
w_hid = random.rand(2,2)*sigw        # [inp1,inp2] x-> [hid1,hid2]
b_hid = random.rand(1,2)*sigw        # 1.0 -> [b_hid1,b_hid2]
w_out = random.rand(2,1)*sigw        # [hid1,hid2] x-> [out1]
b_out = random.mrand(1,1)*sigw        # 1.0 -> [b_out1]
w_out_prev_change = np.zeros(shape(w_out))
b_out_prev_change = np.zeros(shape(b_out))
w_hid_prev_change = np.zeros(shape(w_hid))
b_hid_prev_change = np.zeros(shape(b_hid))

maxepochs = 1000
errors = np.zeros((maxepochs,1))
N = 0.01 # learning rate parameter
M = 0.10 # momentum parameter
