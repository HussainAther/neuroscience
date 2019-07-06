function plotSpikeRasterAndHistogram(spikeTimesMS, plotSynchTimeMS, preSynchDurationMS, postSynchDurationMS, widMS)
%
%  make a raster plot from a set of spike times and a histogram of all the
%  runs (i.e. repeats) that are included:
%
%  This function assumes that the spike times from each run are in a separate cell array:
%  i.e:
%
%  spikeTimesMS{1} = [ 3 10  30 34 ...];    % the list of spike tims on run 1
%  spikeTimesMS{2} = [ 5 9  22 70 ...];     % the list of spike tims on run 2
%   ...
%
%  plotSynchTimeMS = the time in the movie that will be called 0 on the
%       plot (e.g. the start of a new condition)
%  preSynchDurationMS = the amount of time to show before time 0 (should be 0 or greater) 
%  postSynchDurationMS = the amount of time to show AFTER time 0 (should be 0 or greater) 
%  widMS = the wid of each of the bins in the histogram -- try varying this!
%
%

figure; clf;
subplot('position',[0.2 0.15 0.7 0.35]); hold on;

nRuns = length(spikeTimesMS);
axis([-preSynchDurationMS,postSynchDurationMS,0 nRuns+1]);
allSpikeTimesInrange = [];

for run = 1:nRuns
    synchedSpikeTims = spikeTimesMS{run} - plotSynchTimeMS;
    z = find((synchedSpikeTims>=-preSynchDurationMS) & (synchedSpikeTims<=postSynchDurationMS));
    for spikeTime = synchedSpikeTims(z);
		plotTick(run,spikeTime);
    end
    allSpikeTimesInrange = [allSpikeTimesInrange synchedSpikeTims(z)]; 
end
plot([0 0],[0 nRuns+1],'b-');
xlabel('Time from synch (ms)');
ylabel('Run number');


%% make histogram
centersPost = [(widMS/2):widMS:postSynchDurationMS+(widMS/2)];
centersPre = [-(widMS/2):(-widMS):-preSynchDurationMS-(widMS/2)];
centers = [fliplr(centersPre) centersPost];
[counts,x] = hist(allSpikeTimesInrange,centers);
subplot('position',[0.2 0.55 0.7 0.4]); hold on;
r = ((counts/nRuns)/widMS)*1000;        %% spikes/sec
h = bar(centers,r,1);
set(h,'FaceColor',[0.5 0.5 0.5]);
axis([-preSynchDurationMS,postSynchDurationMS,0 (max(r)*1.1)]);
ylabel('spikes/sec');
set(gca,'XTick',[]);    % no x numbers

return;


function plotTick(y,t)
    h = plot([t t],[y-0.4 y+0.4],'k-');
    %set(h,'linewidth',4)
return
    
    
    



