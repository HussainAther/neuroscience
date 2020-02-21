% Simulate three Ornstein-Uhlenbeck processes
%

%relaxation time constant
tau = sqrt(2); %ms

%diffusion parameter
c = sqrt(2); %mu S

sig = sqrt(tau*c/2); % equal to 1;

%time step and time vector
dt = 0.05; %ms
t_v = 0:dt:100;
n_t_v = length(t_v);

%SD relaxation
std_v = (1-exp(-2*t_v/tau));

%numerical simulation factors
ef_1 = exp(-dt/tau);
ef_2 = sqrt(1-exp(-2*dt/tau));

%
%first simulation
%

%initial value
x0 = 0;

%random white noise increment
w_v1 = normrnd(0,1,1,n_t_v);
u_v1 = sig*ef_2*w_v1;

%save space for random vector
x_v1 = zeros(1,n_t_v);

%initial condition
x_v1(1) = x0;

for i =2:n_t_v
    %stochastic update
    x_v1(i) = x_v1(i-1)*ef_1 + u_v1(i);
end;

h_f1 = figure;
h_a1 = subplot(3,1,1);
line('Parent',h_a1,'XData',t_v,'YData',x_v1);
line('Parent',h_a1,'XData',t_v,'YData',std_v,'LineStyle',':');
line('Parent',h_a1,'XData',t_v,'YData',-std_v,'LineStyle',':');
title(' \tau = c = \surd 2','interpreter','default','Parent',h_a1);
set(h_a1,'YLim',[-4 4]);

%
%second simulation
%
%relaxation time constant
tau = sqrt(2)/3; %ms

%diffusion parameter
c = 3*sqrt(2); %mu S

%the ratio c/tau = 9
sig = sqrt(tau*c/2); % equal to 1;

%SD relaxation
std_v = (1-exp(-2*t_v/tau));

%numerical simulation factors
ef_1 = exp(-dt/tau);
ef_2 = sqrt(1-exp(-2*dt/tau));

%initial value
x0 = 0;

%random white noise increment
w_v1 = normrnd(0,1,1,n_t_v);
u_v1 = sig*ef_2*w_v1;

%save space for random vector
x_v1 = zeros(1,n_t_v);

%initial condition
x_v1(1) = x0;

for i =2:n_t_v
    %stochastic update
    x_v1(i) = x_v1(i-1)*ef_1 + u_v1(i);
end;

h_a2 = subplot(3,1,2);
line('Parent',h_a2,'XData',t_v,'YData',x_v1);
line('Parent',h_a2,'XData',t_v,'YData',std_v,'LineStyle',':');
line('Parent',h_a2,'XData',t_v,'YData',-std_v,'LineStyle',':');
title(' \tau = \surd 2/3, c = 3 \surd 2','interpreter','default','Parent',h_a2);
set(h_a2,'YLim',[-4 4]);

%
%third simulation
%
%relaxation time constant
tau = 3*sqrt(2); %ms

%diffusion parameter
c = sqrt(2)/3; %mu S

%the ratio c/tau = 1/9
sig = sqrt(tau*c/2); % equal to 1;

%SD relaxation
std_v = (1-exp(-2*t_v/tau));

%numerical simulation factors
ef_1 = exp(-dt/tau);
ef_2 = sqrt(1-exp(-2*dt/tau));

%initial value
x0 = 0;

%random white noise increment
w_v1 = normrnd(0,1,1,n_t_v);
u_v1 = sig*ef_2*w_v1;

%save space for random vector
x_v1 = zeros(1,n_t_v);

%initial condition
x_v1(1) = x0;

for i =2:n_t_v
    %stochastic update
    x_v1(i) = x_v1(i-1)*ef_1 + u_v1(i);
end;

h_a3 = subplot(3,1,3);
line('Parent',h_a3,'XData',t_v,'YData',x_v1);
line('Parent',h_a3,'XData',t_v,'YData',std_v,'LineStyle',':');
line('Parent',h_a3,'XData',t_v,'YData',-std_v,'LineStyle',':');
title(' \tau = 3 \surd 2, c = \surd 2/3','interpreter','default','Parent',h_a3);
set(h_a3,'YLim',[-4 4]);

%print(handles.figure1,'-depsc2','ou_f1.eps');
