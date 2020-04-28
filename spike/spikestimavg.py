# Spike-triggered stimulus average

spikes = mat["rho"]
stimVals = mat["stim"]
n_timeSteps = 150
spikeTimes = spikes[n_timeSteps:].nonzero()[0] + n_timeSteps
