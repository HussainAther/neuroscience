

"""
When we observe spike train of a sensory neuron, we learn about many different aspects of the stimulus.
This way, we gain information. Shannon entropy is closely related to information. We can think about information transmission by sensory neurons in terms of
various parts of the spike train as part of stimulus space.


"""

def totalProb(variableProbs):
    """
    The additivity of information mean that for any distribution of N independent variables,
    the probability of each of them occuring is the product of each of the individual variables.
    """
    return product(variableProbs)
