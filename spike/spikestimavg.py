# Spike-triggered stimulus average

spikes = mat["rho"]
stimVals = mat["stim"]
n_timeSteps = 150
spikeTimes = spikes[n_timeSteps:].nonzero()[0] + n_timeSteps
sepInterval = np.array([2,20,40,60,80,100])
sepInterval = (sepInterval/2).astype("int")
print("Time intervals:", sepInterval, "steps")

