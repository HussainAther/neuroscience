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

def corrCoef(i, j):
    """
    Calculate the len(i)x(len(j) matrix of pairwise Pearson's correlation coefficients between all combinations of N binned
    spike trains. For each pair of spike trains (equal-sized arrays i and j), we determine the correlation by binning i and j
    at the desired bin size. Let b_i and b_j denote the binary vectors and m_i and m_j their respective averages.
    """
    m_i = sum(i)/len(i) # average of a
    m_j = sum(j)/len(j) # average of b
    b_i = np.where(i > 0.5, 1, 0) # binary vector of i
    b_j = np.where(j > 0.5, 1, 0) # binary vector of i
    result = np.zeros((1,len(i),len(i)))
    for a in range(len(i)):
        for b in range(len(j)):
            num = np.multiply(b_i[a]-m_i[a], b_j[b]-m_j[b]) # numerator of Pearson correlation coefficient function
            den = np.sqrt(np.multiply(b_i[a]-m_i[a], b_i[b]-m_i[b]) * np.multiply(b_j[a]-m_j[a], b_j[b]-m_j[b])) # denominator
            result[a][b] = num / den
    return result

def covariance(i, j):
    """
    Calculate the len(i)x(len(j) matrix of pairwise covariances between all combinations of N binned spike trains.
    See the comment for the corrCoef function.
    """
    m_i = sum(i)/len(i) # average of a
    m_j = sum(j)/len(j) # average of b
    b_i = np.where(i > 0.5, 1, 0) # binary vector of i
    b_j = np.where(j > 0.5, 1, 0) # binary vector of i
    result = np.zeros((1,len(i),len(i)))
    for a in range(len(i)):
        for b in range(len(j)):
            num = np.multiply(b_i[a]-m_i[a], b_j[b]-m_j[b]) # numerator of covariance
            den = (len(i)) # denominator
            result[a][b] = num / den
    return result

def spikeTTc(i ,j):
    """
    For spike trains i and j, compute the spike time tiling coefficient.
    """
