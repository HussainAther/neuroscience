import elephant.conversion as conv
import elephant.kernels as kernels
import math
import neo
import numpy as np
import quantities as pq
import scipy as sp
import scipy.stats
import scipy.signal
import warnings

from bartlett import bartlett # from bartlett.py import the Bartlett"s method function
from neo import io
from neo.core import SpikeTrain

"""
We can measure the information rate of spike firing using the following method. First, we estimate
the signal from spikes (using the noise-reduction approaches), then we compute noise in the estimated
signal (measuring the random errors and removing them), after that we compute the signal-to-noise-ratio of the
estimated signal, and finally calcualte lower bound to information rate from the signal-to-noise ratio.
"""

class AsciiSpikeTrainIO(BaseIO):
    """
    Class for reading/writing SpikeTrains in a text file.
    Each Spiketrain is a line.
    """
    is_readable = True
    is_writable = True
    supported_objects = [Segment, SpikeTrain]
    readable_objects = [Segment]
    writeable_objects = [Segment]
    has_header = False
    is_streameable = False
    read_params = { Segment: [ ("delimiter", {"value": "\t", "possible": ["\t", " ", ",", ";"]}), ("t_start", {"value": 0., }),]}
    write_params = { Segment: [("delimiter", {"value": "\t", "possible": ["\t", " ", ",", ";"]}),]}
    name = None
    extensions = ["txt"]

    mode = "file"

    def __init__(self, filename=None):
        """
        This class read/write SpikeTrains in a text file.
        Each row is a spiketrain. filename is the filename to read/write
        """
        BaseIO.__init__(self)
        self.filename = filename

    def read_segment(self, lazy=False, delimiter="\t", t_start=0. * pq.s, unit=pq.s, ):
        """
        delimiter is the columns delimiter in file  "\t" or one space or two space or "," or ";"
        t_start is the time start of all spiketrain 0 by default. unit is the unit of spike times,
        can be a str or directly a quantity.
        """
        assert not lazy, "Do not support lazy"
        unit = pq.Quantity(1, unit)
        seg = Segment(file_origin=os.path.basename(self.filename))
        f = open(self.filename, "Ur")
        for i, line in enumerate(f):
            alldata = line[:-1].split(delimiter)
            if alldata[-1] == "":
                alldata = alldata[:-1]
            if alldata[0] == "":
                alldata = alldata[1:]
            spike_times = np.array(alldata).astype("f")
            t_stop = spike_times.max() * unit
            sptr = SpikeTrain(spike_times * unit, t_start=t_start, t_stop=t_stop)
            sptr.annotate(channel_index=i)
            seg.spiketrains.append(sptr)
        f.close()
        seg.create_many_to_one_relationship()
        return seg

    def write_segment(self, segment, delimiter="\t",):
        """
        Write SpikeTrain of a Segment in a txt file. Each row is a spiketrain. segment is the segment to write.
        Only analog signals will be written. delimiter  is the columns delimiter in file  "\t" or one space
        or two space or "," or ";". information of t_start is lost.
        """
        f = open(self.filename, "w")
        for s, sptr in enumerate(segment.spiketrains):
            for ts in sptr:
                f.write("%f%s" % (ts, delimiter))
            f.write("\n")
        f.close()

# Information provided to us by spike-train
def fanofactor(spiketrains):
    """
    Evaluates the empirical Fano factor F of the spike counts of
    a list of `neo.core.SpikeTrain` objects. Given the vector v containing the observed
    spike counts (one per spike train) in the time window [t0, t1], F is defined as:
    F = var(v)/mean(v). The Fano factor is typically computed for spike trains representing the
    activity of the same neuron over different trials. The higher F, the larger
    the cross-trial non-stationarity. In theory for a time-stationary Poisson
    process, F=1. spiketrains is a list of neo.SpikeTrain objects, quantity arrays, numpy arrays or lists.
    It"s spike trains for which to compute the Fano factor of spike counts. Returns fano (float or nan),
    the Fano factor of the spike counts of the input spike trains. If an empty list is specified, or if all spike trains are empty, F:=nan.
    """
    # Build array of spike counts (one per spike train)
    spike_counts = np.array([len(t) for t in spiketrains])

    # Compute FF
    if all([count == 0 for count in spike_counts]):
        fano = np.nan
    else:
        fano = spike_counts.var() / spike_counts.mean()
    return fano

def instantaneous_rate(spiketrain, sampling_period, kernel="auto", cutoff=5.0, t_start=None, t_stop=None, trim=False):
    """
    Estimate instantaneous firing rate by kernel convolution.
    spiketrain is the (neo.SpikeTrain or list of neo.SpikeTrain objects) Neo object
    that contains spike times, the unit of the time stamps and t_start and t_stop of the spike train.
    sampling_period (Time Quantity) is the time stamp resolution of the spike times.
    The same resolution will be assumed for the kernel.
    kernel is the (string) "auto" or callable object of :class:`Kernel` from module
    "kernels.py". Currently implemented kernel forms are rectangular, triangular,
    epanechnikovlike, gaussian, laplacian, exponential, and alpha function. Example: kernel =
    kernels.RectangularKernel(sigma=10*ms, invert=False) The kernel is used for convolution
    with the spike train and its standard deviation determines the time resolution of the instantaneous
    rate estimation. Default: "auto". In this case, the optimized kernel width for the rate estimation
    is calculated according to [1] and with this width a gaussian kernel is constructed. Automatized
    calculation of the kernel width is not available for other than gaussian kernel shapes.
    cutoff (float) determines the cutoff of the probability distribution of the kernel, i.e.,
    the considered width of the kernel in terms of multiples of the standard deviation sigma. Default: 5.0
    t_start (Time Quantity (optional)) is the start time of the interval used to compute the firing rate. If None,
    it"s assumed equal to spiketrain.t_start. Default: None. t_stop (Time Quantity (optional)) is the
    End time of the interval used to compute the firing rate (included). If None assumed equal to spiketrain.t_stop. Default: None
    trim (bool) if False, the output of the Fast Fourier Transformation being a longer vector than the input
    vector by the size of the kernel is reduced back to the original size of the considered time
    interval of the spiketrain using the median of the kernel. if True, only the region of the convolved signal is returned, where
    there is complete overlap between kernel and spike train. This is achieved by reducing the length of the output of the Fast Fourier
    Transformation by a total of two times the size of the kernel, and
    t_start and t_stop are adjusted. Default: False Return rate (neo.AnalogSignal) that contains the rate estimation in unit
    hertz (Hz). Has a property "rate.times" which contains the time axis of the rate estimate. The unit of this property
    is the same as the resolution that is given via the argument "sampling_period" to the function.
    Raise TypeError if `spiketrain` is not an instance of :class:`SpikeTrain` of Neo, if `sampling_period` is not a time quantity,
    if `kernel` is neither instance of :class:`Kernel` or string "auto", if `cutoff` is neither float nor int,
    if `t_start` and `t_stop` are neither None nor a time quantity, if `trim` is not bool.
    Raise ValueError if `sampling_period` is smaller than zero.
    Example:
    kernel = kernels.AlphaKernel(sigma = 0.05*s, invert = True)
    rate = instantaneous_rate(spiketrain, sampling_period = 2*ms, kernel)
    Reference:
    ..[1] H. Shimazaki, S. Shinomoto, J Comput Neurosci (2010) 29:171–182.
    """
    # Merge spike trains if list of spike trains given:
    if isinstance(spiketrain, list):
        _check_consistency_of_spiketrainlist(
            spiketrain, t_start=t_start, t_stop=t_stop)
        if t_start is None:
            t_start = spiketrain[0].t_start
        if t_stop is None:
            t_stop = spiketrain[0].t_stop
        spikes = np.concatenate([st.magnitude for st in spiketrain])
        merged_spiketrain = SpikeTrain(np.sort(spikes), units=spiketrain[0].units,
                                       t_start=t_start, t_stop=t_stop)
        return instantaneous_rate(merged_spiketrain, sampling_period=sampling_period,
                                  kernel=kernel, cutoff=cutoff, t_start=t_start,
                                  t_stop=t_stop, trim=trim)

spiketrains = io.AsciiSpikeTrainIO( filename = "../DATA_001.txt")
arrivaltimes = instantaneous_rate(spiketrain, sampling_period) # Elephant function of the first kind of real order and complex argument
T = 30 # limit of time we integrate over

I = integrate.quad(lambda x: arrivaltimes, 0, T) * integate.quad(fanofactor(spiketrains))

"""
As a strategy for measuring information rate of stimulus signals, we can generally
1. Estimate signal from spikes
2. Compute noise in estimate
3. Compure signal to noise ratio of estimates
4. Calculate lower bound to information rate from signal to noise ratio.
"""

def nfs(t, eps=1e-9):
    """
    Noise from spike. We can remove the highest spike from an array t of spikes.
    Eps (epsilon) is accuracy.
    """
    n = len(t)
    x = np.arange(n)
    c = np.argmax(t) # find the max of our spike train array
    d = abs(np.diff(t)) # find the differences (intervals) between the spikes
    try:
        l = c - 1 - np.where(d[c-1::-1]<eps)[0][0] # find the locations in that we include in our accuracy from the left side of the max
        r = c + np.where(d[c:]<eps)[0][0] + 1 # on the right side
    except:
        return False # no spike found
    if (r-l) <= 3: # if the gap around the spike is too small to make meaningful intervals
        l -= 1
        r += 1
    s = int(round((r-l)/2)) # interval around which we make gaps 
    lx = l - s
    rx = r + s
    xgap = np.concatenate((x[lx:l], x[r:rx])) # make gaps at the spike areas
    tgap = np.concatenate((t[lx:l], t[r:rx]))
    z = np.polyfit(xgap, tgap, x) # quadratic fit the gapped array 
    p = np.poly1d(z) # get the p polynomial fit in one dimension
    t[l:r] = p(x[l:r]) # convert our array t accordingly. 
    return t 
  
t = [1,1,1,1,1,10,1,1,1,1,15] # sample spike train
signal = nfs(t) # 1. Estimate signal from spikes

"""
Other methods of separating the signal from spike trains (raw input data) involve
finding the entropy of the stimulus distribution and the corresponding conditional distribution
averaged over the distribution of spike trains. The entropy is determined by 
the setup of the experiment.
"""

# the power spectral density of a spike train during spontaneous activity is the noise spectrum of the neuron

fr, pd = bartlett(signal) # 2. Compute noise in estimate

snr = signaltonoise(t) # 3. Compure signal to noise ratio of estimates

# Shannon-Hartley theorem gives us a maximum rate at which information can be transmitted over a 
# communications channel of a specified bandwidth in the presence of noise.
# We treat signal as bandwidth in this case.

C = t * math.log2(1 + signal/fr) # maximum data rate C

# Calculate lower bound to information rate from signal to noise ratio.

R_info = (1/2) * np.polyint(math.log2(1 + snr)/(2*np.pi)) # np.polyint is the indefinite integral (antiderivative)  
