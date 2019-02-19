import numpy as np
import scipy as sp

"""
We can measure the information rate of spike firing using the following method. First, we estimate
the signal from spikes (using the noise-reduction approaches), then we compute noise in the estimated
signal (measuring the random errors and removing them), after that we compute the signal-to-noise-ratio of the
estimated signal, and finally calcualte lower bound to information rate from the signal-to-noise ratio.
"""

# Information provided to us by spike-train
