import numpy as np

"""
Stochastic models of single ion channels. Transition rate between two chemical states has dimensions 
s^-1. If a ligand must be bound for an ion channel to open, we have three states (as Castillo and
Katz described CK mechanism): two shut states and one open state. State 1 is open with an agonist
molecule bound to a receptor, State 2 has a molecule bound but the channel is shut, and State 3 has
the channel shut with receptor unoccupied. Square matrix Q denotes the transition rates with m x m 
dimensions for m states. Each row must sum to zero. This is a homogenous Markov process. 
"""

# Qf-matrix in which alpha is the rate from 1 to 2, kneg is rate from 2 to 3, 
# kpos is the rate from 3 to 2, and beta is the rate from 2 to 1. xA is the 
# ligand concentration. 
QCK = np.matrix([[-alpha, alpha, 0], 
                 [beta, -(kneg+beta), kneg],
                 [0, kpos*xA, -kpos*xA]])

# Next, we suppose the open channel in the CK model can be blocked by a molecule of a 
# blocker substance with concentration xB of blocker substance. Now there are two shut
# states and the open state has one more shut state. 
