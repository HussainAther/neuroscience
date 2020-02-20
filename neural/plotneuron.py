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
        
    def leakyIF(self):
        '''
        Leaky Integrate and Fire model.
        '''
        self.timeseries = np.linspace(0, self.T-1, self.T)
        self.V = np.ones(self.T)*self.vR
        ii = 0 # index counter
        while ii < self.T-2:
            dV = -(self.gL*(self.V[ii] - self.vR) + self.currentInj[ii])/self.C
            self.V[ii+1] = self.V[ii] + dV
            if self.V[ii+1] >= self.vT:
                self.V[ii+1] = 20
                self.V[ii+1] = self.vR
            ii += 1
       
    def plotNeuron(self):
        '''
        Plot and save the Integrate and Fire (IF) voltage response.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.plot(self.timeseries, self.currentInj, c='k')
        ax.sef_title('Current injection', style='italic')
        ax.set_ylabel('Current (nA)', style='italic')
        ax2 = fig.add_subplot(212)
        ax2.plot(self.timeseries, self.V, c='k')
        ax2.set_title('Integrate and fire voltage response', style='italic')
        ax2.set_xlabel('Time (ms)', style='italic')
        ax2.set_ylabel('Voltage (mV)', style='italic')
        plt.tight_layout()
        plt.savefig('IFVresponse.png')

myNeuron = Neuron()
myNeuron.createInjectionCurrent()
myNeuron.leakyIF()
myNeuron.plotNeuron()
