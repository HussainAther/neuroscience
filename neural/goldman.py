import numpy as np

from fractions import Fraction

'''
Goldman-Hodgkin-Katz (goldman hodgkin katz) equation.
'''

def ghk(K = 310, KO = 4, KI = 140, NaO = 150, NaI = 15, ClO = 0, ClI = 0, PK = 1, PNa = 0.333333333, PCl = 0):
    RT = 8.3144598 * K
    F = 9.6485*10*10*10*10
    up = PK*KO+PNa*NaO+PCl*ClI
    down = PK*KI+PNa*NaI+PCl*ClO
    res = RT/F * (np.log(up/down))
    return('{} Volt'.format(res))
