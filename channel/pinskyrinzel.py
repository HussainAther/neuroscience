"""
Pinsky-Rinzel (pinsky rinzel) model two-compartment biophysical model. It mimics a full, detailed
model of pyramidal cells.

CmdVs(t)/dt = −ILeak(Vs)−INa(Vs,h) − IK−DR(Vs,n) + g (Vd(t)−Vs(t))/p

CmdVd(t)/dt = −ILeak(Vd)−ICa(Vd,s)−IK−AHP(Vd,q)-IK-CVs(t)−Vd(t) − I/(1-p) +gc (Vs(t)-Vd(t))/(1-p)

[Ca]' = -.002ICa - .0125[Ca]
"""

def pr(Ca, Vd,  
