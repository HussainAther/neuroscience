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

# plot the network's performance over the course of learning
# inputs will be a regular grid of [inp1,inp2] points
n_grid = 20
g_grid = np.linspace(-1.0, 2.0, n_grid)
g1,g2 = np.meshgrid(g_grid, g_grid)
plt.figure()

# train
for i in range(maxepochs):
    net_out = np.zeros(shape(xor_out))
    for j in range(shape(xor_in)[0]): # for each training example
        # forward pass
        act_inp = xor_in[j,:]
        act_hid = tansig( (act_inp * w_hid) + b_hid )
        act_out = tansig( (act_hid * w_out) + b_out )
        net_out[j,:] = act_out[0,:]

        # error gradients starting at outputs and working backwards
        err_out = (act_out - xor_out[j,:])
        deltas_out = multiply(dtansig(act_out), err_out)
        err_hid = deltas_out * transpose(w_out)
        deltas_hid = multiply(dtansig(act_hid), err_hid)

        # update the weights and bias units
        w_out_change = -2.0 * transpose(act_hid)*deltas_out
        w_out = w_out + (N * w_out_change) + (M * w_out_prev_change)
        w_out_prev_change = w_out_change
        b_out_change = -2.0 * deltas_out
        b_out = b_out + (N * b_out_change) + (M * b_out_prev_change)
        b_out_prev_change = b_out_change

        w_hid_change = -2.0 * transpose(act_inp)*deltas_hid
        w_hid = w_hid + (N * w_hid_change) + (M * w_hid_prev_change)
        w_hid_prev_change = w_hid_change
        b_hid_change = -2.0 * deltas_hid
        b_hid = b_hid + (N * b_hid_change) + (M * b_hid_prev_change)
        b_hid_prev_change = b_hid_change
