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
[t1, u1]=ode45('rnn_ode', tspan, u0, [], nn, tau_inv, dx, beta, alpha, w, I_ext);
u0=u1(size(t1,1),:0; tspan=[10,20];
I_ext=zeros(nn,1);
[t2, u2]=ode45('rnn_ode', tspan, u0, [], nn, tau_inv, dx, beta, alpha, w, I_ext);
u0=u2(size(t2,1),:); tspan=[20,40];
w=(ws-w_inh).*(1+wa(:, :, 2));
[t3, u3]=ode45('rnn_ode', tspan, u0, [], nn, tau_inv, dx, beta, alpha, w, I_ext);
u0=u3(size(t3,1),:); tspan=[40,50];
w=ws-w_inh; I_ext=zeros(nn,1);
[t4, u4]=ode45('rnn_ode', tspan, u0, [], nn, tau_inv, dx, beta, alpha, w, I_ext);
u0=u4(size(t4,1),:); tspan=[50,70];
w=(ws-w_inh).*(1+2*wa(:,:,1)); I_ext=zeros(nn,1);
[t5, u5]=ode45('rnn_ode', tspan, u0, [], nn, tau_inv, dx, beta, alpha, w, I_ext);
u0=u5(size(t5,1),:); tspan=[70,80];
w=ws-w_inh; I_ext=zeros(nn,1);
[t6, u6]=ode45('rnn_ode', tspan, u0, [], nn, tau_inv, dx, beta, alpha, w, I_ext);

% Plot.
u=[u1; u2; u3; u4; u5; u6]; t=[t1; t2; t3; t4; t5; t6];
r=1./(1+exp(-beta.*(u-alpha)));
%surf(t', 1:nn.r', 'linestyle', 'none')
h=surf(t',1:nn,r');
set(h, 'linestyle', 'none');
