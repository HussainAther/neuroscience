import numpy as np
import matplotlib.pyplot as plt

class Neuron():
    
    def __init__(self):
        self.C = .281  # Capacitance
        self.gL = .030 
        self.vR = -60.6 
        self.vT = -50.4
    
    def createInjectionCurrent(self):
        self.currentInj = np.append(np.zeros(10), np.arange(100)/100)
        self.T = len(self.currentInj)  
