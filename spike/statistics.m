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

t_vect = t_On:dt:(NumTimePoints-1)*dt; % time vector for each trial

% Plot rasters.
figure(1)
for trial=1:NumTrials
   plot(t_vect,trial*spikes(2,:,trial), "+")
   hold on
end
xlabel("time (ms)")
ylabel("trial number")
hold off

% Plot averaged (PSTH) data for one particular angle.
TrialSum_vect = sum(spikes(2,:,:),3); % add up spike trains across trials
SmoothingWidth = 61; % smooth data over +/- [(this # - 1)/2] bins
TrialSum_smooth_vect = smooth(TrialSum_vect,SmoothingWidth);
figure(2)
plot(t_vect,1000*TrialSum_smooth_vect/(NumTrials*dt))
xlabel("time (ms)")
ylabel("Average firing rate (Hz)")
