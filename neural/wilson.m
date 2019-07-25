% Integration of Wilson model with multistep ode solver
clear; clf;

% parameters of model; 1=K,R 2=Ca,T 3=KCa,H 4=Na
g(1)=26; g(2)=2.25; g(3)=9.5; g(4)=1; 
E(1)=-.95; E(2)=1.20; E(3)=E(1); E(4)=.50;
tau(1)=1/4.2; tau(2)=1/14; tau(3)=1/45; tau(4)=1;
tau=tau';
