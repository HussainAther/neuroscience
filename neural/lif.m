% Leaky Integrate-and-fire (IAF) neuron
clear; clf;

% parameters
tau_inv=.1; % inverse time constant
uu=0; % initial membrane voltage
I_ext=12; % constant external input
tspan=[0 100] % integration interval
theta=10; % firing threshold

% Integration with Euler method
t_step=0; for it=0:100;
    t_step=t_step+1;
    x=uu<theta;
    uu=x*(1-tau_inv)*uu+tau_inv*Iext;%+randn;
    u(t_step)=uu;
    t(t_step)=it;
    s(t_step)=1-x;
end

subplot('position',[0.13 0.13 1-0.26 0.6])
    plot(t,u);
    hold on; plt([0 100], [10 10], '--'');
    axis([0 100 0. 12]);
    xlabel('time [\tau]');
    ylabel('u(t)');

subplot('position',[0.13 0.8 1-0.26 0.1])
    plot(t,s,'.', 'markersize', 20);
    axis([0 100 0.5 1.5]);
    set(gca, 'xtick', [], 'ytick', []);
    ylabel('spikes');
