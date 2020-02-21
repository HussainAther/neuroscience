% Frequency response of the passive isopotential cell

function ffreq

A = 4*pi*(1e-6);    % (cm)^2
Cm = 1;             % muF/(cm)^2
GCl = 0.3;          % mS/(cm)^2

tau = Cm/GCl;       %in msec
tau_s = tau*1e-3;   %converted to sec
f = 0:0.1:100;      %Hz

%because Cm is in muF/cm2 and A is in cm2, rin ins in MOhms

rin = 1./(A*Cm*sqrt( (2*pi*f).^2 + 1/tau_s^2 ));

h = figure(1);
set(h,'Position',[1 1 712 358]);

h1 = subplot(1,2,1);
plot(f,rin,'color','k');
set(h1,'TickDir','out');
xlabel('\omega (Hz)','fontsize',14)
ylabel('R_i_n (M\Omega)','fontsize',14)
box off

ph = -atan(2*pi*f*tau_s);
ph_deg = ph*180/pi;
h2 = subplot(1,2,2);
plot(f,ph_deg,'color','k');
xlabel('\omega (Hz)','fontsize',14)
ylabel('phase (deg)','fontsize',14)
set(h2,'TickDir','out');
box off

