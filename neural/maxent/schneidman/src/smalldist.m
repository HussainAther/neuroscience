% Find the marginals (firing rates and correlations) of
% sample data, and compare them against the maximum entropy
% model predictions.

% Load example spiking data of 15 neurons.
load example15

% Randomly divide it into a training set and a test set 
% (so we can verify how well we trained).
[ncells,nsamples] = size(spikes15);
idx_train = randperm(nsamples,ceil(nsamples/2));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes15(:,idx_train);
samples_test = spikes15(:,idx_test);

% Create a k-pairwise model (pairwise maxent with synchrony constraints).
model = maxent.createModel(ncells,"kising");
