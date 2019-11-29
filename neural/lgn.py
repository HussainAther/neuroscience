import numpy as np

"""
Lateral geniculate nucleus (LGN) neurons for ocular dominance,
the way cells receive more input fomr one eye. 
"""
eps = 5 # the correlation between activity in the two populations
        # is epsilon (eps) times as strong as the mean squared activity
        # in each population
C = np.array([[1, eps], [eps, 1]])
