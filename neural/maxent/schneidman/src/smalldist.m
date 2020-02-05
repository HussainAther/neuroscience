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

% Train the model to a threshold of 1.5 standard deviations from the 
% error of computing the marginals. Because the distribution is larger 
% (50 dimensions) we cannot explicitly iterate over all 5^20 states
% in memory and will use markov chain monte carlo (MCMC) methods to obtain an approximation
model = maxent.trainModel(model,samples_train,"threshold",1.5);

% Get the marginals (firing rates and correlations) of the test data and see 
% how they compare to the model predictions. Here the model marginals could not be 
% computed exactly so they will be estimated using monte-carlo. We specify the
% number of samples we use so that their estimation will have the same amoutn noise as the empirical marginal values
marginals_data = maxent.getEmpiricalMarginals(samples_test,model);
marginals_model = maxent.getMarginals(model,"nsamples",size(samples_test,2));

% Plot them on a log scale.
figure
loglog(marginals_data,marginals_model,"b*");
hold on;
minval = min([marginals_data(marginals_data>0)]);
plot([minval 1],[minval 1],"-r"); % identity line
xlabel("empirical marginal");
ylabel("predicted marginal");
title(sprintf("marginals in %d cells",ncells));

% The model that the MCMC solver returns is not normalized. If we want to compare the 
% predicted and actual probabilities of individual firing patterns, we will need to 
% first normalize the model. We will use the wang-landau algorithm for this. We chose 
% parameters which are less strict than the default settings so that we will have a faster runtime.
disp("Normalizing model...");
model = maxent.wangLandau(model,"binsize",0.1,"depth",15);

