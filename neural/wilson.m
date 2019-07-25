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

function ydot=wilson(y, t, flag, I_ext, g, E, tau)
% odefile for wilson equations
% parameters of the model; 1=K,R; 2=Ca,T; 3=KCA,H; 4=Na;
    V=y(4); y3=y(1:3);
    x0(1)=1.24 + 3.7*V + 3.2*V^2;
    x0(2)=4.025 + 11.6*V + 8*V^2;
    x0(3)=3*y(2);
    x0=x0';
    ydot=tau.*(x0-y3);
    y3(4)=17.8+47.6*V+33.8*V^2;
    I=g.*y3.*(y(4)-E);
    ydot(4)=I_ext-sum(I);
return
