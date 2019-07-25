% Lorenz system trajectory attractor
clear; clf;
a=10; b=28; c=8/3;
u0=zeros(3,1)+.4; param=0; tspan=[0, 100];
[t, u]=ode45('lorenz_odefile', tspan, u0, [], a, b, c);
plot3(i(:, 1), u(:, 2), u(:, 3))
