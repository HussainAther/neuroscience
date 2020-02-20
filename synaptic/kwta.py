import numpy as np

'''
We use the k-winners-take-all (kWTA) form of inhibition that uses
a simple set of approximations to the net effects of the inhibitory interneurons (in a
sort of 'winner takes all' approach). We capture the 'set point' nature of inhibitory
feedback as analogous to the desired temperature for the unit. In this example, k refers to the number
of neurons that we choose to become active. 
'''

l = range(10) # original list of values
k = 1 # some arbitrary k

theta = 10 # action potential threshold at which a neuron will fire an action
           # potential output to another neuron

j = int(len(l)*k) # number of elements whose maximum we want 
threshold = sorted(l)[j] # the jth element in the sorted list.
result = np.zeros(len(l)) # initialize zero array 
result = [int(x>=threshold) for x in l] # each list member that is above or equal to the threshold is 1. All others remain 0.
