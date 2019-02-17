import numpy as np
import matplotlib.pyplot as plt

"""
Fourier transformation is a weighted summation of sinusoidal waves that reproduces the original oscillating signal
used in collecting neuroimaging data. The transformed frequency spectrum plot is the graph of the amount of each frequency
found. A uniform signal can have a timestep in between each measurements such that the esries of values can represented as x_k
where k is the time index (from 0 to N-1). The Fast fourier transformation determines the discrete fourier transform
when N is a power of 2. This reduces the number of operations from O(N^2) to O(NlogN).
"""

def savePlot(x, y, xlabel, ylabel):
    plt.plot(x, y, color="k")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    fileName = ylabel.replace(" ","")
    plt.savefig(fileName)
    plt.close()

freqs = np.fft(sig)

freqReal = [f.real for f in freqs] # frequency of real numbers
savePlot(times, freqReal, "freq", "FT real")

freqImag = [f.imag for f in freqs] # frequency of imaginary numbers
savePlot(times, freqImag, "freq", "FT imag")

powerSpec = [abs(f)**2 for f in freqs] # power series of the spectrum
savePlot(times, powerSpec, "freq", "FT power")
