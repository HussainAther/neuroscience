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
