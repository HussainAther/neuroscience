import numpy as np

"""
We use the k-winners-take-all (kWTA) form of inhibition that uses
a simple set of approximations to the net effects of the inhibitory interneurons (in a
sort of "winner takes all" approach). We capture the "set point" nature of inhibitory
feedback as analogous to the desired temperature for the  unit.
"""

l = range(10) # original list of values

theta = 10 # action potential threshold at which a neuron will fire an action
        # potential output to another neuron


