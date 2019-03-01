from scipy.integrate import odeint

"""
In neuroscience, we use kinetmatic  transformations to determine sensory-motor questions
of how the nervous system can choose a single configuartion of body parts (such as arms and legs) in space.
We can use mathematical methods of transformation and processing input values themselves to model
kinematics transformations. We can study how the central nervous system if the predictable
patterns of variability in limb movement can be attributed to these mathematical coordinate transformations.

We define the states of a neuron using ordinary differential equations.
"""

# Set up some parameters that we can work with.
E_params = {
    "E_leak" : -7.0e-2,
    "G_leak" : 3.0e-09,
    "C_m"    : 3.0e-11,
    "I_ext"  : 0*1.0e-10
}

Na_params = {
    "Na_E"          : 5.0e-2,
    "Na_G"          : 1.0e-6,
    "k_Na_act"      : 3.0e+0,
    "A_alpha_m_act" : 2.0e+5,
    "B_alpha_m_act" : -4.0e-2,
    "C_alpha_m_act" : 1.0e-3,
    "A_beta_m_act"  : 6.0e+4,
    "B_beta_m_act"  : -4.9e-2,
    "C_beta_m_act"  : 2.0e-2,
    "l_Na_inact"    : 1.0e+0,
    "A_alpha_m_inact" : 8.0e+4,
    "B_alpha_m_inact" : -4.0e-2,
    "C_alpha_m_inact" : 1.0e-3,
    "A_beta_m_inact"  : 4.0e+2,
    "B_beta_m_inact"  : -3.6e-2,
    "C_beta_m_inact"  : 2.0e-3
}

K_params = {
    "k_E"           : -9.0e-2,
    "k_G"           : 2.0e-7,
    "k_K"           : 4.0e+0,
    "A_alpha_m_act" : 2.0e+4,
    "B_alpha_m_act" : -3.1e-2,
    "C_alpha_m_act" : 8.0e-4,
    "A_beta_m_act"  : 5.0e+3,
    "B_beta_m_act"  : -2.8e-2,
    "C_beta_m_act"  : 4.0e-4
}

Ca_params = {
    "E_Ca" : 150e-03 ,
    "Ca_act_alpha_A" : 0.08e+06,
    "Ca_act_alpha_B" : -10e-03,
    "Ca_act_alpha_C" : 11e-03,
    "Ca_act_beta_A"  : 0.001e+06,
    "Ca_act_beta_B"  : -10e-03,
    "Ca_act_beta_C"  : 0.5e-03,
    "G_Ca"           : 1.0e-08, # G_Ca (uS) ***** this appears as zero in Ekeberg 1991 *****
    "Ca_rho"         : 4.0e+03,
    "Ca_delta"       : 30.0,
    "G_KCA"          : 0.01e-06
}

params = {
    "E_params"  : E_params,
    "Na_params" : Na_params,
    "K_params"  : K_params,
    "Ca_params" : Ca_params
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


    # external current (from "voltage clamp", other compartments, other neurons, etc)
    I_ext = Epar["I_ext"]

    # calculate Na rate functions and I_Na. from Ekeberg, et al., 1991.
    # Na activation
    alpha_act = Na["A_alpha_m_act"] * (E-Na["B_alpha_m_act"]) / (1.0 - exp((Na["B_alpha_m_act"]-E) / Na["C_alpha_m_act"]))
    beta_act = Na["A_beta_m_act"] * (Na["B_beta_m_act"]-E) / (1.0 - exp((E-Na["B_beta_m_act"]) / Na["C_beta_m_act"]) )
    dmdt = ( alpha_act * (1.0 - m) ) - ( beta_act * m )
    # Na inactivation
    alpha_inact = Na["A_alpha_m_inact"] * (Na["B_alpha_m_inact"]-E) / (1.0 - exp((E-Na["B_alpha_m_inact"]) / Na["C_alpha_m_inact"]))
    beta_inact  = Na["A_beta_m_inact"] / (1.0 + (exp((Na["B_beta_m_inact"]-E) / Na["C_beta_m_inact"])))
    dhdt = ( alpha_inact*(1.0 - h) ) - ( beta_inact*h )
    
    # Na-current:
    I_Na =(Na["Na_E"]-E) * Na["Na_G"] * (m**Na["k_Na_act"]) * h


    # calculate K rate functions and I_K
    alpha_kal = K["A_alpha_m_act"] * (E-K["B_alpha_m_act"]) / (1.0 - exp((K["B_alpha_m_act"]-E) / K["C_alpha_m_act"]))
    beta_kal = K["A_beta_m_act"] * (K["B_beta_m_act"]-E) / (1.0 - exp((E-K["B_beta_m_act"]) / K["C_beta_m_act"]))
    dndt = ( alpha_kal*(1.0 - n) ) - ( beta_kal*n )
    
    # K current
    I_K = (K["k_E"]-E) * K["k_G"] * n**K["k_K"]


    # Ca rate functions and Ca current
    alpha_Ca_act = (Ca["Ca_act_alpha_A"]*(E-Ca["Ca_act_alpha_B"]))/(1-exp((Ca["Ca_act_alpha_B"]-E)/Ca["Ca_act_alpha_C"]))
    beta_Ca_act = (Ca["Ca_act_beta_A"]*(Ca["Ca_act_beta_B"]-E))/(1-exp((E-Ca["Ca_act_beta_B"])/Ca["Ca_act_beta_C"]))
    dqdt = alpha_Ca_act*(1-q) - beta_Ca_act*q
    
    # Ca current
    I_Ca = (Ca["E_Ca"] - E)*Ca["G_Ca"]*(q**5)


    # Ca2+ gated K channels
    dCaAPdt = (Ca["E_Ca"] - E)*Ca["Ca_rho"]*(q**5) - Ca["Ca_delta"]*CaAP
    E_K = K["k_E"]
    # Ca2+ gated K current
    I_KCA = (K["k_E"] - E)*Ca["G_KCA"]*CaAP


    # leak current
    I_leak = (Epar["E_leak"]-E) * Epar["G_leak"]

    # calculate derivative of E
    dEdt = (I_leak + I_K + I_Na + I_ext + I_Ca + I_KCA) / Epar["C_m"]
    statep = [dEdt, dmdt, dhdt, dndt, dqdt, dCaAPdt]

    return statep
