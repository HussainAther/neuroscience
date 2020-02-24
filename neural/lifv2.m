% Finite difference method
% Leaky linear integrate and fire method
% Sub-threshold regime
% SI units unless specified
% INPUTS -----------------------------------------------------------------
R = 1e4;                 % membrane resistance
C = 1e-8;                % membrane capacitance
I0 = 2e-3;               % external current amplitude
N = 500;                 % number of data points for calculations


% SETUP ------------------------------------------------------------------
tau = R * C;             % time constant
tmin = 0;                % time interval for modelling
tmax = 10*tau;
t = linspace(tmin, tmax, N);
dt = t(2) - t(1);

Iext = ones(1,N);         % initialize membrane current
V = zeros(1,N);           % initialize membrane voltage - numerical
Va = V;                   % initialize membrane voltage - analytical
K = dt / tau;           % constant

% external stimulus current
Iext = I0 .* Iext;      % constant input


% numerical: fintie difference method
for c = 1 : N - 1
   V(c+1) = V(c) - K * (V(c) - R * Iext(c));
end

% analytical solution
Va = (R * I0) .* (1 - exp(-t ./ tau));

figure(1)
fs = 12;
set(gcf,'units','normalized');
set(gcf,'position',[0.02 0.40 0.3 0.3]);
set(gca,'fontsize',fs);

sfx = 1e3; sfy = 1; x = sfx .* t;   y = sfy .* V;
plot(x,y,'linewidth',2);

xlabel('time t (ms)','fontsize',fs);
ylabel('membrane potential  V  (V)','fontsize',fs);
hold on

Np = (1:10:N); sfx = 1e3; sfy = 1; x = sfx .* t(Np);
y = sfy .* V(Np);
plot(x,y,'or');


