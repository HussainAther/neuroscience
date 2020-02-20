% Train a time-dependent model.
% Load spiking data of 15 neurons.
load ../data/example15_spatiotemporal

ncells = size(spikes15_time_dependent,1);
history_length = 2;

% Create joint words of (x_t-2,x_t-1,x_t).
xt = [];
for i = 1:(history_length+1)
   xt = [xt;spikes15_time_dependent(:,i:(end-history_length-1+i))];
end

% Create a spatiotemporal model that works on series of binary words: (x_t-2,x_t-1,x_t).
% This essentially describes the probability distribution as a second-order Markov process.
% We will model the distribution with a composite model that uses firing rates, total synchrony in the last 3 time bins
% and pairwise correlations within the current time bin and pairwise correlations between the current activity of a cell
% and the previous time bin.
time_dependent_ncells = ncells*(history_length+1);
inner_model_indep = maxent.createModel(time_dependent_ncells,'indep');  % firing rates
inner_model_ksync = maxent.createModel(time_dependent_ncells,'ksync');  % total synchrony in the population

% Add pairwise correlations only within the current time bin.
second_order_correlations = num2cell(nchoosek(1:ncells,2),2);
temporal_matrix = reshape(1:time_dependent_ncells,[ncells,history_length+1]);
temporal_interactions = [];

% Add pairwise correlations from the current time bin to the previous time step.
for i = 1:(history_length)
    temporal_interactions = [temporal_interactions;temporal_matrix(:,[i,i+1])];
end

% Bunch of all this together into one probabilistic model.
temporal_interactions = num2cell(temporal_interactions,2);
inner_model_pairwise = maxent.createModel(time_dependent_ncells,'highorder',[second_order_correlations;temporal_interactions]);
mspatiotemporal = maxent.createModel(time_dependent_ncells,'composite',{inner_model_indep,inner_model_ksync,inner_model_pairwise});

% Train the model on the concatenated words
disp('training spatio-temporal model...');
mspatiotemporal = maxent.trainModel(mspatiotemporal,xt);

% Sample from the model by generating each sample according to the history.
% For this we need to generate the 3n-dimensional samples one by one, each time fixing two-thirds of the code word
% corresponding to time t-2 and t-1 and sampling only from time t.
disp('sampling from spatio-temporal model...');
x0 = uint32(xt(:,1));
xspatiotemporal = [];
nsamples = 10000;
for i = 1:nsamples

    % Get next sample. we will use burn-in of 100 each step to ensure that we don't introduce time-dependent stuff
    % related to the sampling process itself.
    xnext = maxent.generateSamples(mspatiotemporal,1,'fix_indices',1:(ncells*history_length),'burnin',100,'x0',x0);
    generated_sample = xnext(((ncells*history_length)+1):end,1);

    % Shift the 'current' state by one time step.
    x0 = [x0((ncells+1):end,:);generated_sample];

    % Add the current output.
    xspatiotemporal = [xspatiotemporal,generated_sample];

end

% Plot the result.
display_begin = 2000;
nsamples_to_display = 300;

% Plot the actual raster.
pos = [400,600,700,250];
figure('Position',pos);
subplot(2,1,1);
pos = get(gca, 'Position');pos(1) = 0.055;pos(3) = 0.9;set(gca, 'Position', pos);
imshow(~spikes15_time_dependent(:,display_begin+(1:nsamples_to_display)));
title('Actual data');

% Plot raster sampled from a spatiotemporal model.
% Output saved to 'doc/timedep.png'.
subplot(2,1,2);
pos = get(gca, 'Position');pos(1) = 0.055;pos(3) = 0.9;set(gca, 'Position', pos);
imshow(~xspatiotemporal(:,display_begin+(1:nsamples_to_display)));
title('Synthetic data (2nd-order Markov)');
