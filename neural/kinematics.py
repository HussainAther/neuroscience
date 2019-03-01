from scipy.integrate import odeint

"""
In neuroscience, we use kinetmatic  transformations to determine sensory-motor questions
of how the nervous system can choose a single configuartion of body parts (such as arms and legs) in space.
We can use mathematical methods of transformation and processing input values themselves to model
kinematics trnasformations. We can study how the central nervous system if the predictable
patterns of variability in limb movement can be attributed to these mathematical coordinate transformations.

We define the states of a neuron using ordinary differential equations.
"""

def neuron(state, t, params):
    """
    Following the properties of Ekeberg, et al., 1991, we can create a model of a neuron
    to simulate neural networks using the soma and the dendritic tree.
    """

    E = state[0]    # soma potential
    m = state[1]    # Na activation
    h = state[2]    # Na inactivation
    n = state[3]    # K activation
    q = state[4]    # Ca activation
    CaAP = state[5] # Ca2+ dependent K channel

