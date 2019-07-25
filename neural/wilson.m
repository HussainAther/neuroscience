% Integration of Wilson model with multistep ode solver
clear; clf;

% parameters of model; 1=K,R 2=Ca,T 3=KCa,H 4=Na
g(1)=26; g(2)=2.25; g(3)=9.5; g(4)=1; g=g';
E(1)=-.95; E(2)=1.20; E(3)=E(1); E(4)=.50;
tau(1)=1/4.2; tau(2)=1/14; tau(3)=1/45; tau(4)=1;
tau=tau';

% Integration:
%1: Equilibration; no external input;
y0=zeros(4,1); y0(4)=-1; param=0; I_ext=0; tspan=[0 100];
[t, y]=ode15s('wilson_odefile', tspan, y0, [], I_ext, g, E, tau);

%2: Integration with external input;
y0=y(size(t,1),:); param=0; I_ext=1; tspan=[0 200];
[t, y]=ode15s('wilson_odefile', tspan, y0, [], I_ext, g, E, tau);

plot(t, 100*y(:, 4));
