import numpy as np
from scipy.io import loadmat

"""
We can use Gaussians in one dimension to model spatial receptive fields
with center-surround structure.
"""
x = loadmat("ganglion.mat") # read in the MATLAB data using scipy's function loadmat.

list(x) # show all column names
len(x["stim"]) # number of spikes in the whole experiment.
x["__header__"] # view info about the MATLAB file.
x["counts"] # view the counts value for each experiment.
max(x["stim"]) # max value of the stimulus.
min(x["stim"]) # min value of the stimulus.
x["stim"][len(x)/2:] # second half of experiments
x["stim"][:len(x)/2] # first half of experiments
np.var(x["stim"]) # variance
