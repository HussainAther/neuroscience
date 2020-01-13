%{
Spike train statistics
%}

load Sur_Orientation_SpikeData
spikes = double(spikes);

% Initialize variables.
NumControls = 2; % number of control experiments with no grating
dt = 1; % spacing between sampled time points [ms]
t_On = 0; % time stimulus turns on [ms]
t_Move = 500; % time stimulus begins moving [ms]
t_Off = 2500; % time stimulus turns off [ms]
NumAngles = size(spikes,1) - NumControls; % number of angles tested, equally spaced
 %last 2 sets of recordings are controls
NumTimePoints = size(spikes,2); % number of time points time was sampled every 1 ms
NumTrials = size(spikes,3)l % number of trials performed at each angle

t_vect = t_On:dt:(NumTimePoints-1)*dt; % time vector for each trial

% Plot rasters.
figure(1);
for trial=1:NumTrials
   plot(t_vect,trial*spikes(2,:,trial), "+")
   hold on
end
xlabel("time (ms)");
ylabel("trial number");
hold off;

% Plot averaged (PSTH) data for one particular angle.
TrialSum_vect = sum(spikes(2,:,:),3); % add up spike trains across trials
SmoothingWidth = 61; % smooth data over +/- [(this # - 1)/2] bins
TrialSum_smooth_vect = smooth(TrialSum_vect,SmoothingWidth);
figure(2);
plot(t_vect,1000*TrialSum_smooth_vect/(NumTrials*dt));
xlabel("time (ms)");
ylabel("Average firing rate (Hz)");

% Compute Fano factor.
TStartCount = 600; % time to start computing average
TEndCount = 2500; % time to end computing average
NumCounts_vect = sum(spikes(ThisOrientation,((TStartCount+1)/dt):(TEndCount/dt),:),2); # number of counts of spikes per trial
FanoFactor = (std(NumCounts_vect)^2)/mean(NumCounts_vect);

% Compute the interspike intervals (ISI).
SpikeTimes_vect = dt*find(abs(spikes(ThisOrientation,TStartCount:TEndCount,1)-1) < 0.00000001);
isi_vect = diff(SpikeTimes_vect);

% Plot ISI histogram.
figure(3);
hist(isi_vect,8);
xlabel("isi (ms)");
ylabel("Number of occurrences");

% Compute coefficient of variation (CV) ISI.
mean_isi = mean(isi_vect);
std_isi = std(isi_vect);
CV_isi = mean_isi/std_isi;

% Compute average rate for one trial.
TotalCounts = sum(NumCounts_vect); % total spikes across time and trials
AveRate = 1000*TotalCounts/((TEndCount-TStartCount)*NumTrials); % average across time and trials

% Compute tuning curves.
TuningCurveCounts_matrix = sum(spikes(:,((TStartCount+1)/dt):(TEndCount/dt),:),3);
TuningCurveCounts_vect = sum(TuningCurveCounts_matrix,2); %sum over time dimension
TuningCurve_AveRate_vect = 1000*TuningCurveCounts_vect/((TEndCount-TStartCount)*NumTrials);

% Plot tuning curve.
DeltaTheta = 360/NumAngles;
Orientation_vect = 0:DeltaTheta:(NumAngles-1)*DeltaTheta;
figure(4);
plot(Orientation_vect,TuningCurve_AveRate_vect(1:NumAngles));
xlabel("orientation (deg)");
ylabel("Ave firing rate (Hz)");
Control_AveRate = mean(TuningCurve_AveRate_vect(17:18));

% Plot difference from the control's average rate.
TuningCurve_RateDiff_vect = TuningCurve_AveRate_vect - Control_AveRate;
figure(5);
plot(Orientation_vect,TuningCurve_RateDiff_vect(1:NumAngles));
xlabel("orientation (deg)");
ylabel("Ave firing rate relative to control(Hz)");
