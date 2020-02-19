import numpy as np

"""
Izhikevich model for firing
"""

class IzhNeuron:
    def __init__(self, label, a, b, c, d, v0, u0=None):
        self.label = label
        self.a = a
        self.b = b
        self.c = c
        self.d = d
       
        self.v = v0
        self.u = u0 if u0 is not None else b*v0
