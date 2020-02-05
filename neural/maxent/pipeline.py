from genSamples import generateSamples

"""
Pipeline for Maximum entropy analysis
"""

# Generate 5000 samples from the distribution with default starting  
# point of 0 and default burn-in of 10000 samples.
samples = generateSamples(model, 5000)

