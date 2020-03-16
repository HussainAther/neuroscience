
"""
Sample from a maximum entropy model.
"""

def generateSamples(model, nsamples, burnin=1000):
    """
    Accept a model and number of samples to be generated via Metropolis-Hastings (M-H MH).
    Drop the first 1000 generated samples (burn in) and skip most of the samples generated along the 
    way. 
    """
    
