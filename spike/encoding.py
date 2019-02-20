import numpy as np
import matplotlib.pyplot as plt

"""
Encoding models address the question of how the firing rate of a neuron is related to the stimulus that is presented.



A biological example would be where we stimulate a neuron with a particular light intensity.
If the firing rate depends linearly on this light intensity then the first property tells us
that doubling (or halving) the light intensity should double (or halve) the firing rate of the neuron.
The second property says that adding two light intensities together should produce a firing rate that
is equal to the sum of the firing rates to each individual lights. This assumption of linearity does
not hold really well in reality as neural firing rates cannot increase indefinitely.
"""
