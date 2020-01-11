%{
Spike train statistics
}%

load Sur_Orientation_SpikeData
spikes = double(spikes);

% Initialize variables.
NumControls = 2; % number of control experiments with no grating
dt = 1; % spacing between sampled time points [ms]
t_On = 0; % time stimulus turns on [ms]
t_Move = 500; % time stimulus begins moving [ms]
t_Off = 2500; % time stimulus turns off [ms]
NumAngles = size(spikes,1) - NumControls % number of angles tested, equally spaced;
 %last 2 sets of recordings are controls
NumTimePoints = size(spikes,2) % number of time points; time was sampled every 1 ms
NumTrials = size(spikes,3) % number of trials performed at each angle
