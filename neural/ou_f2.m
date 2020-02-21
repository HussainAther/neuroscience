% First simulation ornstein uhlenbeck 
%

%mean excitatory conductance
ge0 = 0.012; %mu S 

%relaxation time constant and steady-state SD
tau_e = 2.7; %ms
sig_e = 0.003; %mu S

%time step and time vector
dt = 0.5; %ms
t_v = 0:dt:1000;
n_t_v = length(t_v);

%numerical simulation factors
ef_1 = exp(-dt/tau_e);
ef_2 = sqrt(1-exp(-2*dt/tau_e));

%initial value of the random conductance variable
x0 = 0;

%random white noise increment
w_ve = normrnd(0,1,1,n_t_v);
u_ve = sig_e*ef_2*w_ve;

%save space for random vector
x_ve = zeros(1,n_t_v);

%initial condition
x_ve(1) = x0;

for i =2:n_t_v
    %stochastic update
    x_ve(i) = x_ve(i-1)*ef_1 + u_ve(i);
end;

%add steady-state value
x_ve = x_ve + ge0;

h_f1 = figure;
h_a1 = subplot(3,2,1);
line('Parent',h_a1,'XData',t_v,'YData',x_ve);
set(h_a1,'YLim',[0 0.04]);
xlabel(h_a1,'time (ms)');
ylabel(h_a1,'conductance (\mu S)');

%histogram of values
[ne, xoute] = hist(x_ve,20);

h_a2 = subplot(3,2,3);
bar(h_a2,xoute,ne,'r');

ye = normpdf(xoute,ge0,sig_e);
ye = (sum(ne)/sum(ye)) * ye;

line('Parent',h_a2,'XData',xoute,'YData',ye);

set(h_a2,'TickDir','out','Xlim',[0 0.04]);
xlabel(h_a2,'conductance (\mu S)');

%power spectrum
f = 0:0.1:200;
pfe = 2*sig_e^2*tau_e*1e-3./(1 + (2*pi*f*tau_e*1e-3).^2);
h_a3 = subplot(3,2,5);
line('Parent',h_a3,'XData',f,'YData',pfe);
xlabel(h_a3,'frequency (Hz)');
ylabel(h_a3,'Power density (\mu S^2 / Hz)');

%%
%second simulation
%

%mean inhibitory conductance
gi0 = 0.057; %mu S 

%relaxation time constant and steady-state SD
tau_i = 10.5; %ms
sig_i = 0.0066; %mu S

%time step and time vector
dt = 0.5; %ms
t_v = 0:dt:1000;
n_t_v = length(t_v);

%numerical simulation factors
if_1 = exp(-dt/tau_i);
if_2 = sqrt(1-exp(-2*dt/tau_i));

%initial value of the random conductance variable
x0 = 0;

%random white noise increment
w_vi = normrnd(0,1,1,n_t_v);
u_vi = sig_i*if_2*w_vi;

%save space for random vector
x_vi = zeros(1,n_t_v);

%initial condition
x_vi(1) = x0;

for i =2:n_t_v
    %stochastic update
    x_vi(i) = x_vi(i-1)*if_1 + u_vi(i);
end;

%add steady-state value
x_vi = x_vi + gi0;

h_a4 = subplot(3,2,2);
line('Parent',h_a4,'XData',t_v,'YData',x_vi);
set(h_a4,'YLim',[0.04 0.08]);
xlabel(h_a4,'time (ms)');
ylabel(h_a4,'conductance (\mu S)');

h_a5 = subplot(3,2,4);
%histogram of values
[ni, xouti] = hist(x_vi,20);
bar(h_a5,xouti,ni,'r');

yi = normpdf(xouti,gi0,sig_i);
yi = (sum(ni)/sum(yi)) * yi;

line('Parent',h_a5,'XData',xouti,'YData',yi);

set(h_a5,'TickDir','out','Xlim',[0.02 0.1]);
xlabel(h_a5,'conductance (\mu S)');

%power spectrum
f = 0:0.1:200; %frequency in Hz
pfi = 2*sig_i^2*tau_i*1e-3./(1 + (2*pi*f*tau_i*1e-3).^2);
h_a6 = subplot(3,2,6);
line('Parent',h_a6,'XData',f,'YData',pfi);
xlabel(h_a6,'frequency (Hz)');
ylabel(h_a6,'Power density (\mu S^2 / Hz)');

%print(handles.figure1,'-depsc2','ou_f2.eps');


