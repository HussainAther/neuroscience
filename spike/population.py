"""
We study the population activity by taking into account the spikes of a population of neurons 
that are sent off to another population (from m to n). The proportion of active neurons in the
presynaptic population m can be used to define population activity.
"""

def delta(n):
    """
    Dirac delta function.
    """
    if n == 0:
        return 1
    else:
        return 0

def activity(t, N):
    """
    Population activity of population of N neurons (array of spikes with 0s and 1s) at time t.  
    Spike train is the sum of Dirac delta functions (using sigma).
    """
    deltat = .5 # some time interval
    sigma = 0
    for i in N:
        sigma += delta(i)
    return sigma / (deltat*len(N)) 
