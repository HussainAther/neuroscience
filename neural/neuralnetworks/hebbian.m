function net = hebbRNN_create_model(N, B, I, p, g, dt, tau, varargin)

% net = hebbRNN_create_model(N, B, I, p, g, dt, tau, varargin)
%
% This function initializes a recurrent neural network for later training
% and execution
%
% INPUTS:
%
% N -- the number of recurrent neurons in network
%
% B -- the number of outputs
%
% I -- the number of inputs
%
% p -- the sparseness of the J (connectivity) matrix, (range: 0-1)
%
% g -- the spectral scaling of J
%
% dt - the integration time constant
%
% tau - the time constant of each neuron
%
%
% OPTIONAL INPUTS:
%
% actFun -- the activation function used to tranform activations into
% firing rates
% Default: 'tanh'
%
% netNoiseSigma - the variance of random gaussian noise added at each time
% point
% Default: 0
%
% useBiasNeurons -- whether or not to included neurons that have a fixed
% value
% Default: false
%
% numBiasNeurons -- how many bias neurons to use
% Default: 0
%
% feedback -- whether or not to feed the output of the network back
% Default: false
%
% energyCost -- how much to weight the firing rate of the output units in
% the calculation of total error. This parameter encourages the network to
% find solutions with as low firing rate as possible.
% Default: 0
%
% outputUnitIdent -- user-defined indexing of output units
% Default: 1:B
%
% biasUnitIdent -- user-defined indexing of bias units
% Default: []
%
%
% OUTPUTS:
%
% net -- the network structure

actFunType = 'tanh'; % Default activation function
netNoiseSigma = 0.0; % Default noise-level
useBiasNeurons = false; % Default use of bias neurons
numBiasNeurons = 0; % Default number of bias neurons
feedback = false; % Default use of output feedback
energyCost = 0; % Default weight of cost function
outputUnitIdent = 1:B; % Default identity of output units
biasUnitIdent = []; % Default identity of bias units
co(1:2) = false; % Both unit identities must be provided if defaults are not used
optargin = size(varargin,2);

for i = 1:2:optargin
    switch varargin{i}
        case 'actFun'
            actFunType = varargin{i+1};
        case 'netNoiseSigma'
            netNoiseSigma = varargin{i+1};
        case 'useBiasNeurons'
            useBiasNeurons = varargin{i+1};
        case 'numBiasNeurons'
            numBiasNeurons = varargin{i+1};
        case 'feedback'
            feedback = varargin{i+1};
        case 'energyCost'
            energyCost = varargin{i+1};
            
        case 'outputUnitIdent'
            outputUnitIdent = varargin{i+1};
            co(1) = true;
        case 'biasUnitIdent'
            biasUnitIdent = varargin{i+1};
            co(2) = true;
    end
end

%% Assertions
assert(co(1) == co(2), 'Both identities must be provided.')
assert(islogical(feedback), 'Must be logical.')
assert(p >= 0 && p <= 1, 'Sparsity must be between 0 and 1.')

%% Set bias unit indices
if useBiasNeurons
    if isempty(biasUnitIdent)
        net.biasUnitIdent = B+1 : B+numBiasNeurons;
    else
        net.biasUnitIdent = biasUnitIdent;
    end
else
    net.biasUnitIdent = [];
end
