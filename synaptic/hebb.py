import numpy as np

"""

"When the axon of cell A is near enough to excite a cell B and repeatedly or persistently
takes part in firing it, some growth process or metabolic change takes place in one or
both cells such the A's efficiency, as one of the cells firing B, is increased." - Psychologist Donald Hebb in "The Organization of Behavior"

Long term increases in potentiation in synaptic strength can be called the LTP (long term potentiation).

"""
def deltaW(r_i, q_j):
    """
    Change in weight connecting neuron j to neuron i.
    r_i and q_j are the firing rates of neurons i and j.
    alpha (learning rate) determines how much the weight changes during each association.
    phi_pre is the presynpatic LTD threshold.
    phi is the postsynaptic LTD threshold.
    """
    alpha = 5
    phi_pre = 10
    phi = 12
    return alpha * (q_j - phi_pre) * (r_i - phi)


"""
We can express the rule using vectors and find an outer product. In this case,
r and q are vector components. Sensory-motor matching can explain learning associations
between different sensory modalities. We can sum the outer product of the vectors r and
q to determine the vector of the weight matrix.
"""

def deltaWvec(r, q):
    """
    Returns a weight matrix using the restriction of requireing orthogonal memories
    for learning as the readout of the learning is done using the original coordinate
    system.
    """
    result = 0
    for i in np.outer(r, q):
        result += i
    return result


