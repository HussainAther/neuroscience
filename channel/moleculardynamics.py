import numpy as np

"""
We use the atomistic simulation technique of ion channels as a general example of molecular dyanmics.
We use empirical potentisl to simulate interactions between all atoms in the system with an example
potential function.
"""

def V():
    """
    Example molecular dynamics potential V.
    """
    bondsum = 0 # summation of bond energies
    anglesum = 0 # summation of angle energies
