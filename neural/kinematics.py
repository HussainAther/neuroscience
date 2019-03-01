from scipy.integrate import odeint

"""
In neuroscience, we use kinetmatic  transformations to determine sensory-motor questions
of how the nervous system can choose a single configuartion of body parts (such as arms and legs) in space.
We can use mathematical methods of transformation and processing input values themselves to model
kinematics trnasformations. We can study how the central nervous system if the predictable
patterns of variability in limb movement can be attributed to these mathematical coordinate transformations.

We define the states of a neuron using ordinary differential equations.
"""

# Set up some parameters that we can work with.
E_params = {
    'E_leak' : -7.0e-2,
    'G_leak' : 3.0e-09,
    'C_m'    : 3.0e-11,
    'I_ext'  : 0*1.0e-10
}

Na_params = {
    'Na_E'          : 5.0e-2,
    'Na_G'          : 1.0e-6,
    'k_Na_act'      : 3.0e+0,
    'A_alpha_m_act' : 2.0e+5,
    'B_alpha_m_act' : -4.0e-2,
    'C_alpha_m_act' : 1.0e-3,
    'A_beta_m_act'  : 6.0e+4,
    'B_beta_m_act'  : -4.9e-2,
    'C_beta_m_act'  : 2.0e-2,
    'l_Na_inact'    : 1.0e+0,
    'A_alpha_m_inact' : 8.0e+4,
    'B_alpha_m_inact' : -4.0e-2,
    'C_alpha_m_inact' : 1.0e-3,
    'A_beta_m_inact'  : 4.0e+2,
    'B_beta_m_inact'  : -3.6e-2,
    'C_beta_m_inact'  : 2.0e-3
}

K_params = {
    'k_E'           : -9.0e-2,
    'k_G'           : 2.0e-7,
    'k_K'           : 4.0e+0,
    'A_alpha_m_act' : 2.0e+4,
    'B_alpha_m_act' : -3.1e-2,
    'C_alpha_m_act' : 8.0e-4,
    'A_beta_m_act'  : 5.0e+3,
    'B_beta_m_act'  : -2.8e-2,
    'C_beta_m_act'  : 4.0e-4
}

Ca_params = {
    'E_Ca' : 150e-03 ,
    'Ca_act_alpha_A' : 0.08e+06,
    'Ca_act_alpha_B' : -10e-03,
    'Ca_act_alpha_C' : 11e-03,
    'Ca_act_beta_A'  : 0.001e+06,
    'Ca_act_beta_B'  : -10e-03,
    'Ca_act_beta_C'  : 0.5e-03,
    'G_Ca'           : 1.0e-08, # G_Ca (uS) ***** this appears as zero in Ekeberg 1991 *****
    'Ca_rho'         : 4.0e+03,
    'Ca_delta'       : 30.0,
    'G_KCA'          : 0.01e-06
}

params = {
    'E_params'  : E_params,
    'Na_params' : Na_params,
    'K_params'  : K_params,
    'Ca_params' : Ca_params
}

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

    Epar = params["E_params"]
    Na   = params["Na_params"]
    K    = params["K_params"]
    Ca   = params["Ca_params"]
