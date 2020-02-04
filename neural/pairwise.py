"""
Fit a pairwise modeel to data. Similar to Schneidman, et al. (2006). "Weak pairwise correlations imply strongly correlated network states in a neural population." Nature.

This uses gradient descent on negative log-likelihood with gradient estimated
by averaging over samples from the model.
"""
