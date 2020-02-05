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

% Train the model to a threshold of 1.5 standard deviations from the 
% error of computing the marginals. Because the distribution is larger (50 dimensions) 
% we cannot explicitly iterate over all 5^20 states in memory and will use Markov Chain Monte 
% Carlo (MCMC) methods to obtain an approximation.
model = maxent.trainModel(model,samples_train,"threshold",1.5);

% Get the marginals (firing rates and correlations) of the test data and see how they compare to the model predictions.
% here the model marginals could not be computed exactly so they will be estimated using monte-carlo. We specify the
% number of samples we use so that their estimation will have the same amoutn noise as the empirical marginal values
marginals_data = maxent.getEmpiricalMarginals(samples_test,model);
marginals_model = maxent.getMarginals(model,"nsamples",size(samples_test,2));

% Plot them on a log scale.
% Output saved to "doc/largedist.png".
figure
loglog(marginals_data,marginals_model,"b*");
hold on;
minval = min([marginals_data(marginals_data>0)]);
plot([minval 1],[minval 1],"-r"); % identity line
xlabel("empirical marginal");
ylabel("predicted marginal");
title(sprintf("marginals in %d cells",ncells));
