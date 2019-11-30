import numpy as np

"""
Long short-term memory (long short term lstm) model
of a recurrent neural network.
"""

ct = [0, 0, 0] # candidate layer
ht = [0, 0, 0] # hidden layer
