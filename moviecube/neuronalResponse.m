function [responseIPS] = neuronalResponse(startTimeMS,endTimeMS,spikeTimesMS)
%
%  Return response rate from list of spike times and interval to analyze
%  units = spikes/sec (impulses per sec = IPS)
%
%
if (startTimeMS>= endTimeMS) error('startTimeMS should be greater than endTimeMS'); return; end;
z = find((spikeTimesMS>=startTimeMS) & (spikeTimesMS<endTimeMS));
r = length(z)/(endTimeMS-startTimeMS);      %% spikes per ms
responseIPS = r*1000;                       %% spikes per sec
