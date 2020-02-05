% Use randomly projected models to check the similarity between them
% to determine how well the patterns of activity emerge among first,
% second, and third-order correlations.

% Load spiking data of 15 neurons.
load data/example15

% Randomly divide it into a training set and a test set (so we can verify how well we trained).
[ncells,nsamples] = size(spikes15);
idx_train = randperm(nsamples,ceil(nsamples/2));
idx_test = setdiff(1:nsamples,idx_train);
samples_train = spikes15(:,idx_train);
samples_test = spikes15(:,idx_test);

% Create a random projection model with default settings.
model = maxent.createModel(ncells,"rp");

% Train the model to a threshold of one standard deviation from the error of computing the marginals.
% Because the distribution is relatively small (15 dimensions) we can explicitly represent all 2^15 states
% in memory and train the model in an exhaustive fashion.
model = maxent.trainModel(model,samples_train,"threshold",1);

% Now check the kullback-leibler divergence between the model predictions and the pattern counts in the test-set.
empirical_distribution = maxent.getEmpiricalModel(samples_test);
model_logprobs = maxent.getLogProbability(model,empirical_distribution.words);
test_dkl = maxent.dkl(empirical_distribution.logprobs,model_logprobs);
fprintf("Kullback-Leibler divergence from test set: %f\n",test_dkl);

% Create a random projection model with a specified number of projections and specified sparsity.
model = maxent.createModel(ncells,"rp","nprojections",500,"sparsity",0.1);

% Train the model.
model = maxent.trainModel(model,samples_train,"threshold",1);

% Now check the Kullback-Leibler (KL) divergence between the model predictions and the pattern counts in the test-set.
empirical_distribution = maxent.getEmpiricalModel(samples_test);
model_logprobs = maxent.getLogProbability(model,empirical_distribution.words);
test_dkl = maxent.dkl(empirical_distribution.logprobs,model_logprobs);
fprintf("Kullback-Leibler divergence from test set: %f\n",test_dkl);
