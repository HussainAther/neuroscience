% Integration of Wilson model with the Euler method.
clear; clf;
% parameters of model; 1=K,R 2=Ca,T 3=KCa,H 4=Na
g(1)=26; g(2)=2.25; g(3)=9.5; g(4)=1; 
E(1)=-.95; E(2)=1.20; E(3)=E(1); E(4)=.50;
dt=.01; I_ext=0; V=-1; x=zeros(1,4); tau(1)=dt./4.2;
tau(2)=dt./14; tau(3)=dt./45; tau(4)=1;

% Integration
t_rec=0;

for t=-100:dt:200
    switch t;
        case 0; I_ext=1;
    end
    x0(1)=1.24 + 3.7*V + 3.2*V^2;
    x0(2)=4.205 + 11.6*V + 8*V^2;
    x0(3)=3*x(2);
    x0(4)=17.8 + 47.6*V + 33.8*V^2;
    x=x-tau.*(x-x0); %rem x(4)=x0(4) because tau(4)=1
    I=g.*x.*(V-E);
    V=V+dt*(I_ext-sum(I));
    if t>=0;
        t_rec=t_rec+1;
        x_plot(t_rec)=t;
        y_plot(t_rec)=V;
    end
end % time loop

plot(x_plot, 100*y_plot)
