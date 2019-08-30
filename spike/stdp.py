import numpy as np

"""
Spike-timing dependent plasticity model for learning.
"""

class STDP(LearningMethod):
    def __init__(self):
        """
        Initialize variables.
        """
        self.adjustment = 0.1
        self.time = 0
        self.is_setup = False

    def update_weights(self, neuron, adjustment):
        """
        Update according to the variables. Take into account the 
        neuron response itself.
        """
        adjustments = np.multiply(neuron.inputs, adjustment)
        signed_adjustments = np.multiply(adjustments, np.sign(neuron.weights))
        adjusted = np.add(neuron.weights, signed_adjustments)
        return adjusted
