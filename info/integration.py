import numpy as np
import matplotlib.pyplot as plt

"""
Humans use Bayesian optimal integration in the information from both
visual and auditory stimuli at the same tiem and same location in time and space.

Following "Bayesian inference with probabilistic population codes" (2006) by Ma et al.,
we can express the mathematical distribution of subject's estimates of these stimuli.

In a cue combination task, the goal is to integrate two cues, c1 and c2, both of which provide
information about the same stimulus, s. For instance, s could be the spatial location of a stimulus,
c1 could be a visual cue for the location, and c2 could be an auditory cue.

The mean and variance of posterior mu and sigma combine the auditory and visual inputs.

mu3 = (sigma2^2 / (sigma1^2 + sigma2^2)) * mu1 + (sigma1^2/(sigma1^2+sigma2^2) * mu2

1/(sigma3^2) = 1/sigma1^2 + 1/sigma2^2

"""

def visDist(sigma1, K):
    """
    Compute gain1 of the vision distribution.
    """
    return 1/((sigma1**2) * K)


def auditDist(sigma2, K):
    """
    Compute gain2 of the auditory distribution.
    """
    return 1/((sigma2**2) * K)

def visuAuditDist(sigma3, K):
    """
    Compute gain3 of the visual-auditory distribution.
    """
    return visDist(sigma1, K) + auditDist(sigma2, K)


"""
In "Bayesian inference with probabilistic population codes", Wei Ji Ma and colleagues have
proposed that populations of activity automatically probability distributions and that,
due to Poisson noise, multisensory integration could be realised very simply just by summing
the activity of the two populations of activity.

We describe populations with N = 50 neurons and tuning curves that explain the mean spike count
of each neuron in 1 second as a function of the stimulus direction theta.
"""

f_curve = [] # tuning curve function

for G in range( ):
    
