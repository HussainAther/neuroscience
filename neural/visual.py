import numpy as np

from scipy.io import loadmat
from scipy.stats import norm

np.random.seed(0)

'''
We can use Gaussians in one dimension to model spatial receptive fields
with center-surround structure. We use built-in python libraries to show the Gaussian
integrals.
'''

x = loadmat('ganglion.mat') # read in the MATLAB data using scipy's function loadmat.

list(x) # show all column names
len(x['stim']) # number of spikes in the whole experiment.
x['__header__'] # view info about the MATLAB file.
x['counts'] # view the counts value for each experiment.
max(x['stim']) # max value of the stimulus.
min(x['stim']) # min value of the stimulus.
x['stim'][len(x)/2:] # second half of experiments
x['stim'][:len(x)/2] # first half of experiments
np.var(x['stim']) # variance

# Random Gaussians
gauss_1 = norm(loc=-5,scale=5)
gauss_2 = norm(loc=8,scale=3)
gauss_3 = norm(loc=1.5,scale=1)
