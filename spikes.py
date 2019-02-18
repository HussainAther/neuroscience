import sys, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import patsy as pt
import pymc3 as pm

"""
Simple, approximate descriptions of spike statistics. They can be useful for analytic statements about the
qimplications of whole classes of such models for idfferent interesting quantities such as the reliability of code
or the possibility of decoding the spike train.

We can begin with a Poisson model. One spike firing with some probability per unit time can depend on time
but no explicitly on the occurrence times of the other spikes.
"""

# Poisson theta values
theta_spike_bin = .7 # a spike occurs
theta_nospike = .3 # there is no spike
theta_

# samples
q = 1000
df = pd.DataFrame({)
