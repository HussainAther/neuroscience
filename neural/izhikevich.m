%{
Integrate and Fire (IAF) Izhikevich neurons.
}%
a = .03;
b = .25;
c = -60;
d = 4;
tau = .25;
tspan = 0:tau:1000;
T1 = zeros(size(tspan)); 
T1(dsearchn(tspan',200):dsearchn(tspan',800)) = 1;
% membrane potential
V = V + tau*(.04*V^2 + 5*V + 140 - u + T1(ti));
u = u + tau*a*(b*V-u);
if V > 30% there was a spike
    VV(ti+1)=30;
    V=c; u=u+d;
else % there was no spike
    VV(ti+1)=V;
end
uu(ti+1)=u;
