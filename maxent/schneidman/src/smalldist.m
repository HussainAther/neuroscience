% Find the marginals (firing rates and correlations) of
% sample data, and compare them against the maximum entropy
% model predictions.

% Load example spiking data of 15 neurons.
load ../data/example15

% Randomly divide it into a training set and a test set 
% (so we can verify how well we trained).
[ncells,nsamples] = size(spikes15);
idx_train = randperm(nsamples,ceil(nsamples/2));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes15(:,idx_train);
samples_test = spikes15(:,idx_test);

% Create a k-pairwise model (pairwise maxent with synchrony constraints).
model = maxent.createModel(ncells,'kising');

% Train the model to a threshold of 1.5 standard deviations from the 
% error of computing the marginals. Because the distribution is larger 
% (50 dimensions) we cannot explicitly iterate over all 5^20 states
% in memory and will use markov chain monte carlo (MCMC) methods to obtain an approximation
model = maxent.trainModel(model,samples_train,'threshold',1.5);

% Now check the kullback-leibler divergence between the model predictions and 
% the pattern counts in the test-set.
empirical_distribution = maxent.getEmpiricalModel(samples_test);
model_logprobs = maxent.getLogProbability(model,empirical_distribution.words);
test_dkl = maxent.dkl(empirical_distribution.logprobs,model_logprobs);
fprintf('Kullback-Leibler divergence from test set: %f\n',test_dkl);

model_entropy = maxent.getEntropy(model);
fprintf('Model entropy: %.03f   empirical dataset entropy: %.03f\n', model_entropy, empirical_distribution.entropy);

% Get the marginals (firing rates and correlations) of the test data and 
% see how they compare to the model predictions.
marginals_data = maxent.getEmpiricalMarginals(samples_test,model);
marginals_model = maxent.getMarginals(model);

% Plot them on a log scale.
% Output saved to 'doc/smalldist.png'.
figure
loglog(marginals_data,marginals_model,'b*');
hold on;
minval = min([marginals_data(marginals_data>0)]);
plot([minval 1],[minval 1],'-r'); % identity line
xlabel('empirical marginal');
ylabel('predicted marginal');
title(sprintf('marginals in %d cells',ncells));
