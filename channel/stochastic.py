import numpy as np

"""
Stochastic models of single ion channels. Transition rate between two chemical states has dimensions 
s^-1. If a ligand must be bound for an ion channel to open, we have three states (as Castillo and
Katz described CK mechanism): two shut states and one open state. State 1 is open with an agonist
molecule bound to a receptor, State 2 has a molecule bound but the channel is shut, and State 3 has
the channel shut with receptor unoccupied. 
"""
