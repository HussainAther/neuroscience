import numpy as np

"""
We use the atomistic simulation technique of ion channels as a general example of molecular dyanmics.
We use empirical potentisl to simulate interactions between all atoms in the system with an example
potential function.
"""

def chunks(l, n):
    """
    Return successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]

def potential(k, bonds, angles, torsions, eps, sigma, qi, qj, r):
    """
    Example molecular dynamics potential V. We use harmonic terms for bonds and angles, a cosine expansion
    for torsion angles, and Lennard-Jones and Coulomb interactions for non-bonded interactions. k is an array of the
    harmonic force constants. bonds is the list of bond lengths, angles is the list of angles (in degree), torsions 
    are 3-piece tuple to describe the dihedral angles (rotations around a central bond: n, V, gamma), eps and sigma are Lennard-
    Jones parameters in 2-D array form, qi and qj are partial atomic charges (1-D arrays), and r is the distance between 
    the two atoms (2-D array). 
    """
    bondsum = 0 # summation of bond energies
    anglesum = 0 # summation of angle energies
    torsionsum = 0 # summation of torsion energy due to torsion
    li0 = 1 # reference bond length
    thetai0 = 90 # reference bond angle
    m = len(bonds) # number of inputs
    omega = chunks(range(90), m) # get evenly spaced list of 0 to 90 for the number of inputs
    for i in range(m): # for each input
        V = torsions[i][0] # n for dihedral angle
        n = torsions[i][1] # omega for dihedral angle
        gamma = torsions[i][2] # gamma for dihedral angle
        bondsum += (bonds[i] - li0)**2 *(k[i]/2)
        anglesum += (angles[i] - thetai0)**2 *(k[i]/2)
        torsionsum += (1 + np.cos(n*omega-gamma)) * (V/2i)
    ljcsum = 0 # Lennard-Jones Coulomb sum
    for i in range(m):
        for j in range(i+1, m):
            ljcsum += 4*eps[ij] *(((sigma[i][j]/r[i][j]) **12 - (sigma[i][j]/r[i][j])**6)) + qi[i]qj[j]/(4.np.pi*np.finfo(float).eps*r[i][j])
    return bondsum + anglesum + torsionsum + ljcsum 

"""
We can use the Nernst-Planck equation to describe flux of ions driven by an electrochemical
potential gradient across the ion channel. The flux is:

J = -D_i(r) [deltan_i(r,t) + (n_i(r, t)/kT) deltamu_i(r)]

D_i is diffusion coefficient of species i, n_i is position depdendent number density, q_i is the charge, phi is 
the electrostatic potential.

We can combine it with the Poisson equation to derive:

delta[eps(r)deltaph(r) = -4*pi[rho(r) + summation from i=1 to N of z_i*en(r)]

in which the first right-side term is charge density of the fixed charges, and the second term is average charge density of mobile charges.
"""

def chempot(i):
    """
    Return the chemical potential for for each species i.
    """

def pnp_ss():
    """
    When deltaJ equals zero, we get the steady-state equation for drift-diffusino to accomodate the fluxes of mobile ions. 
    """
    N = len(
    for i in range(N)
