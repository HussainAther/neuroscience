% 1-d continuous attractor neural network for path integration
clear; clf; hold on; view(0, 90);
nn = 100; dx=100/nn; % nodes and resolution in deg
tau_inv = 1./1; % inverse of membrane time constant
beta = .06; alpha=.0; % transfer function 1/(1+exp(-beta*(u-theta))

% weight matriecs
sig = 5/dx;
[ws, wa]=hebb_trace(nn, sig);
w_inh=7*(sqrt(2*pi)*sig)^2/nn;

% external input to initiate bubble, no idiothetic cue, w sym
u0=zeros(nn,1)-10; tspan=[0, 10]; param=0;
w=ws-w_inh;
I_ext=zeros(nn,1); for i=-20:20; I_ext(nn/2+i)=50; end
