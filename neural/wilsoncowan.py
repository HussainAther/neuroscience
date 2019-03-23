import numpy as np

"""
The Wilson-Cowan model describes dynamics of interactions between simple exchitatory
and inhibitory model neurons. It's used in modeling neuronal popuatlions with phase
plane methods and numerical solutions to describe the responses of neuronal populations
to stimuli. H.R. Wilson and Jack D. Cowan created the model in the 1970s to prove the existence
of limit cycle dynamics in response to one class of stimuli implies the existence
of multiple stable states and hysteresis in response to a different class of stimuli.
"""

def sig(x)
    """
    Sigmoid function of some value x.
    """
    return 1 / (1 + np.exp(-x))

def fire(ie0, ie1, ii0, ii1, w, t, dt, uu, vv, wee, wei, wie, wii, ze, zi):
    """
    Each variable represents a different neuron and function in our firing model.
    ie = current of the excitatory neuron
    ii = current of the inhibitory neuron
    dt = time step size
    wee, wie, wii, wei = synaptic weights for excitatory and inhibitory neurons, respsectively.
        They connect the synaptic coupling strength between excitatory and inhibitory
        neurons based on the two letters.
        Use these to create a spatially homogenous weighing function.
    u and v = proportion of exchitatory and inhibitory cells firing, respectively.
    z_e and z_e = constant modulatory currents applied ot hte populations.
    """
    i_e = ie0 + ie1 * np.sin(w*t) # current through excitatory neuron
    i_i = ii0 + ii1 * np.sin(w*t) # current through inhibitory neuron
    dE = dt * (-uu + sig((wee * uu) - (wie * vv) - ze + i_e))/tau # excitatory differential
    dI = dt * (-vv + sig((wei * uu) - (wii * vv) - zi + i_i))/tau # inhibitory differential
    uu_p = uu + dE # convert to scaled number
    vv_p = vv + dI
    return uu_p, vv_p

# initial conditions
tau = 1 # timescale
ie0 = 0 # initial current of the excitatory neuron
ii0 = 0 # initial current of the inhibitory neuron
