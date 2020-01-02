import brian2 as b2
import matplotlib.pyplot as plt
import numpy as np

from neurodynex.hodgkin_huxley import HH
from neurodynex.tools import input_factory

%{
Hodgkin-Huxley (hodgkin huxley) model of squid axon.
}%

HH.getting_started()
