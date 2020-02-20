%{
The autocovariance of Y is C_YY = (2alpha^2/pi) * arcsin(C_XX/sigma_X^2+l^2)
Use this result and Exercise 7 to derive numerically the coherence between X and Y when X is white noise with a cut-off frequency of 100 Hz. Show that this allows you to reproduce the red curve in Figure 18.1D. 
%}

% Coherence estimation
%relaxation time constant and steady-state SD
tau_ou = 2.7; %ms
sig_ou = 0.55; %mu S

%time step and time vector
dt = 0.5; %ms
fs = 1/(dt*1e-3); %sampling freq
fn = fs/2; %nyquist freq

%this gives 16 segments to average over for the 
%power spectrum
n_seg = 32;
seg_len = 2048;
n_v = 0:(n_seg*seg_len-1);
t_v = n_v*dt;
n_t_v = length(t_v)

%numerical simulation factors
nsf_1 = exp(-dt/tau_ou);
nsf_2 = sqrt(1-exp(-2*dt/tau_ou));

%initial value of the random conductance variable
x0 = 0;;
