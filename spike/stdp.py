import numpy as np

"""
Spike-timing dependent plasticity model for learning.
"""
class Neuron:
  def __init__(self, weights=0):
    np.seterr(all="ignore")
    self.input = 0
    self.value = 0
    self.output = 0
    self.threshold = 0
    self.fired = False
    self.potential = 0
    self.weights = np.array([self.init_weight(weights) for x in range(weights)], dtype='float64')

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

    def update_trace(self, neuron):
        """
        Trace the neuron to figure out the pathways of it.
        """
        trace = [ 0 if x == 0 else neuron.trace[int(i)] + 1 for (i, x) in enumerate(1 - neuron.inputs)]
        return trace
  
    def setup(self, neuron):
        """
        Create the neuron and trace.
        """
        neuron.trace = np.zeros(len(neuron.inputs))
