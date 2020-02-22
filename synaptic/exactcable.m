% Graph the exact solution of the cable driven by
%
%  I(x,t) = -q_1(x)*(exp(-t)-exp(-2*t))/500
%

a = 1e-4;
ell = 0.1;
C_m = 1;		% micro F / cm^2
g_L = 1/15; %0.3;     		% mS / cm^2
R_2 = 0.3; % 0.034;		% k Ohm cm

tau = C_m/g_L;
lambda = sqrt(a/(2*R_2*g_L));

zeta = (1+lambda^2*pi^2/ell^2)/tau;

[x,t] = meshgrid(0:.004:ell,0:.2:15);
v = cos(x*pi/ell)*sqrt(2/ell)./(C_m*2*a*pi);
v = v.*(exp(-t)*(zeta-2)+exp(-2*t)*(1-zeta)+exp(-zeta*t));
v = -v./(zeta-1)./(zeta-2)/500;

mesh(x,t,v)
colormap([0 0 0])
axis tight
xlabel('x  (cm)','fontsize',16)
ylabel('t  (ms)','fontsize',16)
zlabel('v  (mV)','fontsize',16)

