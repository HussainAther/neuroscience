% 1-d continuous attractor neural network for path integration
clear; clf; hold on; view(0, 90);
nn = 100; dx=100/nn; % nodes and resolution in deg
tau_inv = 1./1; % inverse of membrane time constant
beta = .06; alpha=.0; % transfer function 1/(1+exp(-beta*(u-theta))

% weight matriecs
sig = 5/dx;
[ws, wa]=hebb_trace(nn, sig);
w_inh=7*(sqrt(2*pi)*sig)^2/nn;
