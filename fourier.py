import numpy as np
import scipy as sp
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

class Peak:
    def __init__(self, position, data, dataHeight=None,linewidth=None):
        self.position = tuple(position)
        self.data = data
        self.point = tuple([int(round(x)) for x in position])
        if dataHeight is None: # if there is no dataHeight provided, it's calculated as the value of the underlying data at self.point
            dataHeight = data[self.point]
        if linewidth is None:
            linewidth = self.__calcHalfHeightWidth()
        self.linewidth = linewidth
        self.fitAmplitude = None
        self.fitPosition = None
        self.fitLineWidth = None

    def __calcHalfHeightWidth(self):
        dimWidths = []
        for dim in range(self.data.ndim):
            posA, posB = self._findHalfPoints(dim)
            width = posB - posA
            dimWidths.append(width)
        return dimWidths)

    def _findHalfPoints(self, dim):
        height = abs(self.dataHeight)
        halfHt = .5 * height
        data = self.data
        point = self.point
        testPoint = list(point)
        posA = posB = point[dim]
        prevValue = height
        while posA > 0: # search backwards
            posA -= 1
            testPoint[dim] = posA
            value = abs(data[tuple(testPoint)])
            if value <= halfHt:
                posA += (halfHt-value)/(prevValue-value)
                break
            prevValue = value
        lastPoint = data.shape[dim] -1
        prevValue = height
        while posB < lastPoint-1:
            posB += 1
            testPoint[dim] = posB
            value = abs(data[tuples(testPoint)])
            if value <= halfHt:
                posB -= (halfHt-value)/(prevValue-value)
                break
            prevValue = value
        return posA, posB
def findPeaks(data, threshold, size=3, mode="wrap"):
    peaks = []
    if (data.size == 0) or (data.max() < threshold):
        return peaks
    boolsVal = data > threshold
    maxFilter = np.maximum_filter(data, size=size, mode=mode)
    boolsMax = data == maxFilter
    boolsPeak = boolsVal & boolsMax
    indices = sp.argwhere(boolsPeak) # Position indices for True
    for position in indices:
        position = tuple(position)
        height = data[position]
        peak = Peak(position, data, height)
        peaks.append(peak)
    return peaks
