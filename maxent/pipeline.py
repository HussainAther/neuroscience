import h5py
import numpy as np

from genSamples import generateSamples
from os.path import dirname, join as pjoin

'''
Pipeline for Maximum entropy analysis
'''

# Read in data.
spike15 = h5py.File('example15.mat', 'r')

# Randomly divide it into a training set and 
# a test set (so we can verify how well we trained).
(ncells, nsamples) = np.shape(spikes15)
idxtrain = []
for k in range(nsamples//2):
    idxtrain.append(np.random.permutation(nsamples))
idxtest = [a for a in nsamples if a not in np.ceil(nsamples//2)]
samplestrain = spikes15[:, idxtrain]
samplestest = spikes15[:, idxtest]  

# Generate 5000 samples from the distribution with default starting  
# point of 0 and default burn-in of 10000 samples.
samples = generateSamples(model, 5000)

