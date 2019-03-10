from scipy.io import loadmat

"""
We're going to use data from experiments described in R. Wessel, C. Koch, and F. Gabbiani, Coding of
time-varying electric field amplitude modulations in a wave-type electric fish. J Neurophysiol 75:2280-93 (1996). The
weakly electric fish Eigenmannia has a special organ that generates an oscillating electric field with a frequency of
several hundred Hz. It also has an electrosensory organ, with which it is able to sense electric fields.
"""

x = loadmat("fish.mat") # read in the MATLAB data using scipy's function loadmat.

len(x["stim"]) # number of spikes in the whole experiment.
x["__header__"] # view info about the MATLAB file.
x["rho"] # view the rho value for each experiment.
 max(x["stim"]) # max value of the stimulus.
