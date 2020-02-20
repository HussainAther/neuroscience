import numpy as np

'''
Fisher informaiton used to summarize hte likelihood.
'''

def hessian(x):
    '''
    Calculate the hessian matrix with finite differences for input array x.
    '''
    x_grad = np.gradient(x) 
    hess = np.empty((x.ndim, x.ndim) + x.shape, dtype=x.dtype) 
    for k, grad_k in enumerate(x_grad):
        # iterate over dimensions
        # apply gradient again to every component of the first derivative.
        tmp_grad = np.gradient(grad_k) 
        for l, grad_kl in enumerate(tmp_grad):
            hess[k, l, :, :] = grad_kl
    return hess

def E(a): 
    '''
    Expectation value E for array a. Return the average value 
    of repetitions of the experiment it represents.
    '''
    n = len(a)
    prb = 1/n  # prob of each element
    summ = 0
    for i in range(0, n): 
        sumn += (a[i] * prb)  
    return float(summ) 

def fisher(l):
    '''
    For expectation with respect to p(N|theta) E and likelihood values l, compute 
    Fisher information. 
    ''' 
    return -E(hessian(l))
