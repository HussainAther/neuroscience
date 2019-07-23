%{
Integrate and Fire (IAF) Izhikevich neurons.
}%
a = .03;
b = .25;
c = -60;
d = 4;
tau = .25;
tspan = 0:tau:1000;
T1 = zeros(size(tspan)); T1(dsearchn(tspan',200):dsearchn(tspan',800)) = 1;
