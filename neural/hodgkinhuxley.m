% Hodkin-Huxley equations (hodgkin huxley) integration with Euler method
clear; clf;
% maximal conductances (mS/cm^2); 1=K, 2=Na, 3=R
g(1)=36
g(2)=120
g(3)=.3
% battery voltage (mV); 1=n, 2=m, 3=h
E(1)=-12
E(2)=115
E(3)=10.613
