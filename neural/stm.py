import numpy as np

'''
Short-term (short term stm) memory equations.
'''

def f(u):
    '''
    Positive feedback weight.
    '''
    return 5*u

def g(u):
    '''
    Negative feedback weight.
    '''
    return 2*u

def additive(x, i), A, B, C:
    '''
    In the additive model, we add the terms, possibly nonlinear, that
    determine the rate of change of neuronal activites, or potentials, x,
    with a single potential term i. Return the rate of change for that
    particular potential. 
    '''
    pd = -A[i]*x[i] # Passive decay term using a nonlinear signal
              # with a weight value.
    pf = 0 # Positive feedback term summed over the individual 
           # potentials
    for j in range(len(x)):
        pf += f(j)*B[j]*z[j]
    nf = 0 # Negative feedback term
    for j in range(len(x)):
        nf += g(j)*C[j]
    return pd + pf - nf + I[i] 
