% Solve the inverse passive fiber for the alpha func synaptic input 
%
% NOTE: At figure 3 the user will be prompted to click at a point
%
% usage:        trapcabsyninv(cab,stim,pinc)
%
%               cab.rad = cable radius (cm)
%               cab.ell = cable length (cm)
%    	        cab.dx = space step (cm)
%               cab.dt = timestep (ms)
%               stim.t1 = stim start time (ms)
%               stim.tau = time constant of alpha function (ms)
%               stim.Gsyn = amplitude of synaptic conductance (mS)
%               stim.loc = location of current pulse (cm)
%               stim.Tfin = stopping time (ms)
%               pinc = number of time steps between plots
% example:      
%               cab = struct('rad',1e-4,'ell',.1,'dx',1e-3,'dt',0.05)
%               stim = struct('t1',1,'tau',1,'Gsyn',5e-6,'loc',0.03,'Tfin',100)
%               pinc = 4 
%

function trapcabsyninv(cab,stim,pinc)

Cm = 1;		% micro F / cm^2
g_L = 1/15; %0.3; %1/15; %0.3;     		% mS / cm^2
R_2 = 0.3; %0.15; %0.3; % 0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;                % patch length
A = 2*pi*cab.rad*dx;            % patch surface area
x = dx/2:dx:cab.ell-dx/2;       % vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

v = zeros(Nx,1);		% initial conditions

e = ones(Nx,1);
S = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = 1/dx/dx;
S(Nx,Nx) = 1/dx/dx;

tau = Cm/g_L;
lambda = sqrt(cab.rad/(2*R_2*g_L));
A = 2*pi*cab.rad*dx;

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
vsyn = 70;
Iapp = stim.Gsyn*(dt/2)/A/Cm;

B = (dt/2)*(speye(Nx)+lambda^2*S)/tau;
B = speye(Nx) + B;
bloc = (eloc-1)*(Nx+1)+1;
dBe = B(bloc);

t = 0;
if pinc
figure(1)
plot3(x,t*e,v,'k')
hold on
end

c0 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);
t = dt;
c1 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);

r = zeros(Nx,1);
r(eloc) = vsyn*(c0 + c1)';

vend = zeros(2,Nt);

for j=2:Nt,
    
    B(bloc) = dBe + c1;

    v = B\r; 

    vend(:,j) = v([1 end]);

    if mod(j,pinc) == 0
        plot3(x,t*e,v,'k')
    end

    t = j*dt;
    c0 = c1;
    c1 = Iapp.*((t-stim.t1)./stim.tau).*exp(1-(t-stim.t1)./stim.tau).*(t>stim.t1);

    r = 2*v - r;
    r(eloc) = r(eloc) + vsyn*(c0 + c1)';

end

if pinc
xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)
end

hold off

v0 = vend(1,:);
v1 = vend(2,:);

% find the location of the synapse via means of end potentials

n = length(v0);
v0n = v0.*(1+2*randn(1,n)/100);
fv0 = fft(v0n);
v1n = v1.*(1+2*randn(1,n)/100);
fv1 = fft(v1n);

figure(1)
t = linspace(0,stim.Tfin,Nt);
plot(t,v0n,'k')
hold on
plot(t,v1n,'r')
hold off
box off
legend('v_0','v_1')
xlabel('t (ms)','fontsize',14)
ylabel('v (mV)','fontsize',14)

figure(2)
T = stim.Tfin;
omega = [0:floor(n/2) -(ceil(n/2)-1:-1:1)]/T;
plot(1e3*omega(1:ceil(n/2)),abs(fv0(1:ceil(n/2)))/n,'k')
box off
xlim([0 1000])
ylim([0 0.1])
xlabel('\omega  (Hz)','fontsize',14)
ylabel('$|\hat {\rm v}_0|$','interpreter','latex','fontsize',14)

dt = cab.dt;
ell = cab.ell;

m0 = sum(v0n)*dt;
m1 = sum(v1n)*dt;

sigma = cosh((x-ell)/lambda)./cosh(x/lambda);
c = m0/m1;		

figure(3)
plot(x,sigma,'k')
hold on
plot(x,c*ones(size(x)),'r')
xlabel('x (cm)','fontsize',16)
ylabel('\sigma','fontsize',16)
disp('Park the crosshairs at the intersection point')
pause(.1)
grid
[xb,cb] = ginput(1);
hold off

% the fft approach

cutoff = 0.35;   % 12 is good for orig pap
mask = ones(size(omega));
ind = find( abs(omega) > cutoff );
mask(ind) = 0;
omega = mask.*omega;
mu = sqrt(1+i*2*pi*tau*omega)/lambda;
rvsyn = real(ifft(fv0.*mask.*cosh(xb*mu)));
rvsyn_x_left = real(ifft(fv0.*mask.*mu.*sinh(xb*mu)));
rvsyn_x_right = -real(ifft(fv1.*mask.*mu.*sinh((ell-xb)*mu)));
%gsyn_rec = lambda^2*(rvsyn_x_right-rvsyn_x_left)./(rvsyn - vsyn);
gsyn_rec = (2*pi*cab.rad*g_L)*lambda^2*(rvsyn_x_right-rvsyn_x_left)./(rvsyn - vsyn);

figure(4)
%trueg = (dx*2/dt/g_L)*Iapp.*((t-stim.t1)./stim.tau).*exp(1-(t-stim.t1)./stim.tau).*(t>stim.t1);
trueg = stim.Gsyn.*((t-stim.t1)./stim.tau).*exp(1-(t-stim.t1)./stim.tau).*(t>stim.t1);
plot(t,1e6*trueg,'k')
hold on
plot(t,1e6*gsyn_rec,'r')
legend('true','recovered')
xlabel('t   (ms)','fontsize',16)
ylabel('g_{syn}  (nS)','fontsize',16)
xlim([0 12])
box off
hold off

