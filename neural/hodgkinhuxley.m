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
% Initialize variables
I_Ext=0
V=-10
x=zeros(1,3)
x(3)=1
t_rec=0
% Time step for integration
dt=.01
%Integration with Euler Method
for t=-30:dt:50
    if t==10; I_ext=10; end
    if t==40; I_ext=0; end
    % alpha and beta functions used by Hodgkin and Huxley
    alpha(1)=(10-V)/(100*(Exp((10-V)/10)-1));
    alpha(2)=(25-V)/(10*(Exp((25-V)/10)-1));
    alpha(3)=.07*exp(-V/20);
    beta(1)=.125*exp(-V/80);
    beta(2)=4*exp(-V/18);
    beta(3)=1/(exp((30-V)/10)+1);
    tau=1/(alpha+beta);
    x_0=alpha.*tau;
    % leaky integration with Euler method
    x=(1-dt./tau).*x+dt./tau.*x_0;
    gmh(1)=g(1)*x(1)^4;
    gmh(2)=g(2)*x(2)*3*x(3);
    gmh(3)=g(3);
    I=gnmh.*(V-E);
    V=V+dt*(I_ext-sum(I));
    if t>=0:
        t_rec=t_rec+1;
        x_plot(t_rec)=t;
        y_plot(t_rec)=V;
    end
end % time loop
plot(x_plot, y_plot);
xlabel('Time');
ylabel('Voltage'); 
