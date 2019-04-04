import numpy as np

"""
Represent spike train  as a function of time to place an infinitely thin, infinitely tall window
to mark each spike. Use the Dirac delta function centered at t = t_hat:

1. delta(t) = 0 for t != 0
2. integral from a to b of delta(t)dt = 1 for a <= 0 <= b.
3. integral from a to b of f(t)*delta(t) = f(0) for any function f(t) and for a <= 0 <= b.

We solve this using a convolution integral of two functions. These smoothing functions can convert
messy, noisy data into spike trains that can be analyzed.
"""

A = [1,2,3,4,5,6,7,8,9]
b = [1,2,3]

np.convolve(A, b) # numpy's built-in convolution operator

"""
Interspike interval (ISI) is the time between spikes. ISI rate can be used to calcualte the number
of spikes over time. Then we can calculate a peri-stimulus time histogram (PSTH) from many trials of the
same experiment.
"""

def interspike_interval_histogram(self, bins):
    """
    Return ISI in discrete time bins bins.
    """
    intervals = self.interspike_intervals()
    return np.diff([np.count_nonzero(intervals < (t / 1000.)) for t in [0] + list(bins)]) / float(intervals.size)


"""
We can take into account the Fano factor with every counting interval.
"""

def fano_factor(self, counting_intervals):
    """
    Compute the Fano factor sigma^2_n / mean_n with every given interval. It measures the 
    dsipersion of a probability distributino of Fano noise. 
    """
    ls = []
    for i in counting_intervals:
        counts = self.spike_counts(i)
        ls.append(np.var(counts) / np.mean(counts))
    return np.asarray(ls)
