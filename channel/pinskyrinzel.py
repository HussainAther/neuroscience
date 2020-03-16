import numpy as np

"""
Pinsky-Rinzel (pinsky rinzel) model two-compartment biophysical model. It mimics a full, detailed
model of pyramidal cells. It"s related to the two-compartment abstract integrate-and-fire model.

CmdVs(t)/dt = −ILeak(Vs)−INa(Vs,h) − IK−DR(Vs,n) + g (Vd(t)−Vs(t))/p

CmdVd(t)/dt = −ILeak(Vd)−ICa(Vd,s)−IK−AHP(Vd,q)-IK-CVs(t)−Vd(t) − I/(1-p) +gc (Vs(t)-Vd(t))/(1-p)

[Ca]" = -.002ICa - .0125[Ca]
"""

def Ileak(V):
    """
    Leak current for a given potential.
    """
    return np.sin(V) 

def INa(V, h):
    """
    With Hodgkin-Huxley (Hodgkin huxley) parameter h and potential, 
    return sodium current.
    """
    return np.sin(V)*h

def IK(V, n):
    """
    Potassium current for potential and n parameter.
    """
    return np.sin(V)*n

def pr(Ca, Vd, Vs, h, n, t, s, q, c, p) 
    """
    Pinksy-Rinzel model basic form for input Ca concentration, somatic membrane potential Vm,
    membrane potential of dendritic component Vd and various other variables with ratio of 
    membrane area of somatic component to whole cell p.
    """
    ICa = 1 # Calcium current
    Caa = -.002*ICa - .0125*Ca # Adjusted Calcium current
    gc = 1 # Capacitance-dependent conductance 
    return Caa
