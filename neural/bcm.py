import numpy as np

"""
During neural development, sensory stimulation induces long-term changes in the receptive 
field of the neurons that encode the stimuli. We use hte Bienenstock-Cooper-Munro (BCM bcm) 
model to analyze this process computationallly. We model the BCM plasticity model for a neuron
with N synapses. 
"""

def bcm(x, y):
    """
    We create an (N+1)-dimensional system with two equations. 
    """
     
