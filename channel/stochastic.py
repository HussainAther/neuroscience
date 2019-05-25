import numpy as np

"""
Stochastic models of single ion channels. Transition rate between two chemical states has dimensions 
s^-1. If a ligand must be bound for an ion channel to open, we have three states (as Castillo and
Katz described CK mechanism): two shut states and one open state. State 1 is open with an agonist
molecule bound to a receptor, State 2 has a molecule bound but the channel is shut, and State 3 has
the channel shut with receptor unoccupied. Square matrix Q denotes the transition rates with m x m 
dimensions for m states. Each row must sum to zero. This is a homogenous Markov process. 

Qf-matrix in which alpha is the rate from 1 to 2, kneg is rate from 2 to 3, 
kpos is the rate from 3 to 2, and beta is the rate from 2 to 1. xA is the 
ligand concentration.
"""
 
QCK3 = np.matrix([[-alpha, alpha, 0], 
                 [beta, -(kneg+beta), kneg],
                 [0, kpos*xA, -kpos*xA]])

"""
Next, we suppose the open channel in the CK model can be blocked by a molecule of a 
blocker substance with concentration xB of blocker substance. Now there are two shut
states and the open state has one more shut state, giving us a 4th state. kposB is the rate from 1 to 4, and
knegB is the rate from 4 to 1.
"""

QCK4 = np.matrix([[-(alpha+kposB*xB), alpha, 0, kposB*xB],
                    [beta, -(beta+kneg), kneg, 0],
                    [0, kpos*xA, -kpos*xA, 0],
                    [knegB, 0, 0, -knegB]])

"""
Finally we suppose that nicotinic acetylcholine receptor may have one agonist or two molecules bound
to the shut receptor or the open receptor. If the channel could open without a bound agonist, we'd
have one more state, state 5. In this diagram, kpos2 is the rate from 1 to 2, alpha2 is from 2 to 3,
2*kneg2 is from 3 to 4, kneg is from 4 to 5, 2*kpos is from 5 to 4, beta1 is from 4 to 1, kpos2 is from
4 to 3, beta2 is from 3 to 2, and 2*kneg2 is from 2 to 1.  
"""

QCK5 = np.matrix([[(alpha1 + kpos2*xA, kpos2*xA, 0, alpha1, 0],
                  [2*kneg2, -(alpha2+2*kneg2), alpha2, 0, 0],
                  [0, beta2, -(beta2+2*kneg2), 2*kneg2, 0],
                  [beta1, 0, kpos2*xA, -(beta1+kpos2*xA+kneg), kneg],
                  [0, 0, 0, 2*kpos*xA, -2*kpos*xA]])

"""
In the absence of energy, a system moves to thermodynamic equilibrium at the same rate in each direction on average 
for each individual reaction. This principle of microscopic reverseibility means if there's a cycle in the reaction,
there can be no tendency to move round the cycle in one particular direction. 
"""
