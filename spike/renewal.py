"""
We use Markov renewal processes to generalize Markov jump processes. Renewal processes are stochastic
point processes to describe a sequence of events in time. In the narrow sense, they presuppose stationary 
input and are defined by the fact that the state of the system, and the probability of generating the
next event depends only on the "age" (t-t_hat) of the system. 
"""

def survivor(t, t_hat):
    """
    The survivor function defines the probability that the neuron stays quiescent between t_hat and t.
    It's the probability a neuron "survives" from t_hat to t without firing.
    """

