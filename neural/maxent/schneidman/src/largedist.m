% Working with a larger distribution of neurons using the Markov Chain
% Monte Carlo approach.

% Load 50-neuron dataset.
load data/example50

% Randomly divide into train/test sets.
[ncells,nsamples] = size(spikes50);
idx_train = randperm(nsamples,ceil(nsamples/2));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes50(:,idx_train);
samples_test = spikes50(:,idx_test);

% Create a pairwise maximum entropy model.
model = maxent.createModel(50,"pairwise");
