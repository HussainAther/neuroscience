import numpy as np

"""
Fisher informaiton used to summarize hte likelihood.
"""

import numpy as np

def hessian(x):
    """
    Calculate the hessian matrix with finite differences for input array x.
    """
    x_grad = np.gradient(x) 
    hess = np.empty((x.ndim, x.ndim) + x.shape, dtype=x.dtype) 
    for k, grad_k in enumerate(x_grad):
        # iterate over dimensions
        # apply gradient again to every component of the first derivative.
        tmp_grad = np.gradient(grad_k) 
        for l, grad_kl in enumerate(tmp_grad):
            hess[k, l, :, :] = grad_kl
    return hess

def fisher(E, 
