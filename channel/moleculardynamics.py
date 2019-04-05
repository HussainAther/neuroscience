import numpy as np

"""
We use the atomistic simulation technique of ion channels as a general example of molecular dyanmics.
We use empirical potentisl to simulate interactions between all atoms in the system with an example
potential function.
"""

def V(k, bonds, angles, torsions, eps, sigma, eq, qi, rij):
    """
    Example molecular dynamics potential V. We use harmonic terms for bonds and angles, a cosine expansion
    for torsion angles, and Lennard-Jones and Coulomb interactions for non-bonded interactions. k is an array of the
    harmonic force constants. bonds is the list of bond lengths, angles is the list of angles (in degree), torsions 
    are 3-piece tuple to describe the dihedral angles (rotations around a central bond), eps and sigma are Lennard-
    Jones parameters, eq and qi are partial atomic charges, and rij is the distance between the two atoms. 
    """
    bondsum = 0 # summation of bond energies
    anglesum = 0 # summation of angle energies
    li0 = 1 # reference bond length
    thetai0 = 90 # reference bond angle
    n = len(bonds) # number of inputs
    for i in range(n): # for each input
        bondsum += (bonds[i] - li0)**2 *(k[i]/2)
        anglesum += (angles[i] - thetai0)**2 *(k[i]/2)
