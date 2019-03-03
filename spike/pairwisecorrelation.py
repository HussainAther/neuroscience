import numpy as numpy
import matplotlib.pyplot as plt

"""
Quantifying the degree of correlation between neural spike trains is a key part of analyses of experimental data in many systems.
Neural coordination is thought to play a key role in information propagation and processing and also in self-organization of the
neural system during development. For example, correlated activity plays a critical role in forming the retinotopic map.
In the developing retina, waves of correlated spontaneous activity in retinal ganglion cells have been recorded [on multielectrode arrays (MEAs) and by calcium imaging]
in vitro in many species (Wong, 1999), and shown in vivo using calcium imaging in mouse. These waves show both temporal and spatial correlations.
Much work has focused on assessing the role of this activity in the development of the retinotopic map; typically, both the map and various statistics of the activity
are compared between wild-type and mutant genotypes. The results are used to make inferences about which features of the activity are implicated in retinotopic map formation.
There is strong evidence that correlation between neuronal spike times is involved in this process.

We can determine the correlation index between two spike trains A and B as the factor by which the firing rate
of neuron A increases over its mean value if measured within a fixed window of spikes frmo neuron B.
"""

def corrCoef(a, b,):
    """
    For each pair of spike trains (equal-sized arrays a and b), we determine the correlation by binning i and j
    at the desired bin size. Let b_i and b_j denote teh binary vectors adn m_i and m_j their respective averages.
    """
    m_i = sum(a)/len(a) # average of a
    m_j = sum(b)/len(b) # average of b
    for i in len(a):
        np.multiply(b_i)
