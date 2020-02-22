% Solve the passive cable with current injection via the trapezoid rule
%
% usage:        phi = trapcabextra(cab,stim)
%
%               cab.rad = cable radius (cm)
%               cab.ell = cable length (cm)
%    	        cab.dx = space step (cm)
%               cab.dt = timestep (ms)
%
%               stim.on = start of current alpha function (ms)
%               stim.tau = time constant of alpha (ms)
%               stim.max = maximum of current alpha (micro amps)
%               stim.loc = location of current pulse (cm)
%               stim.Tfin = stopping time (ms)
%
% example:      
%           cab = struct('rad',1e-4,'ell',.1,'dx',1e-3,'dt',0.01)
%           stim = struct('on',0,'tau',1,'max',3e-4,'loc',0.06,'Tfin',20)
%

function phi = trapcabextra(cab,stim)

Cm = 1;		% micro F / cm^2
g_L = 0.3; % mS / cm^2
R_2 = 0.1; % k Ohm cm
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

B = (speye(Nx)+lambda^2*S)/tau;

[L,U] = lu(speye(Nx)+B*dt/2);

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
Iapp = stim.max*(dt/2)/A/Cm;

Istim = e1;
ycenters = ((1:Nx)'-1/2)*cab.ell/Nx;    % coordinates of compartment centers

t = 0;
figure(1)
plot3(t*e,x*1e4,v,'k')
hold on
figure(2)
plot3(t*e,x*1e4,v,'k')
hold on

f0 = 0;
t = dt;
ett = (t-stim.on)/stim.tau;  % elapsed scaled time
f1 = Iapp.*(ett>0)*ett*exp(1-ett); % alpha stimulus

r = zeros(Nx,1);
r(eloc) = (f0 + f1)';
phi = zeros(4,9,Nt);   % initialize extracellular potential

sig = 0.003; % S/cm, extracellular conductivity, page 14 Magneto ref
V0 = zeros(Nx,1);  % zeroth order moment of v
V1 = zeros(Nx,1);  % first order moment of v

for j=1:Nt,
    
    vpre = v;
    v = U \ ( L \ r );   
    
%     if mod(j,24) == 0
%        plot3(t*e,x,v,'k')
%     end
    
    V0 = V0 + v*dt;     % accumulate moments
    V1 = V1 + t*v*dt;

    t = (j+1)*dt;
    f0 = f1;
    ett = (t-stim.on)/stim.tau;  % elapsed scaled time
    f1 = Iapp.*(ett>0)*ett*exp(1-ett); % alpha stimulus

    r = 2*v - r;
    r(eloc) = r(eloc) + (f0 + f1)';

    % collect the membrane currents

    Istim(eloc) = f1;
    Imem = A*(g_L*v - (2/dt)*Istim + Cm*(v-vpre)/dt); % total membrane current
    if mod(j,24) == 0
       figure(1)
       plot3(t*e,x*1e4,v,'k')
       figure(2)
       plot3(t*e,x*1e4,Imem/A,'k')
    end
    for ix = 1:4
        xrec = 0.005*ix;    % horizontal distance between cable and shank
        for iy = 1:9
            yrec = 0.01*iy; % vertical height of electrode
            dscale = 4*pi*sig*sqrt(xrec^2+(yrec-ycenters).^2);
            phi(ix,iy,j) = sum(Imem./dscale);
        end
    end

end

figure(1)
ylabel('y (\mum)','fontsize',14)
xlabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)
text(2,.09,7,'(A)','fontsize',20)
hold off
figure(2)
ylabel('y (\mum)','fontsize',14)
xlabel('t (ms)','fontsize',14)
zlabel('I_{mem} (\mu A)','fontsize',14)
zlim([-310 5])
%text(2,.09,7,'(B)','fontsize',20)
hold off

tim = (1:Nt)*dt;
figure(3)
for ix = 1:4,
for iy = 1:9,
    plot(5+tim+(5/4)*(ix-1)*Nt*dt,squeeze(phi(ix,iy,:))+0.5*(iy-1),'k','linewidth',1.5)
    hold on
end
end

X = [-2 -1 -1 -2];
Y = [-.5 -.5 4.5 4.5];
fill(X,Y,'k')   % draw cable
for j=1:9
    text(-18,(j-1)/2,[num2str(j*100) ' \mum']); % vertical labels
end
for j=1:4,
    text(5+(5/4)*(j-1)*stim.Tfin,-.3,[num2str(50*j) ' \mum']) % horizontal labels
end

plot([10 15],[3.9 3.9],'k')
text(10,3.75,'5 ms')

plot([11 11],[2.1 2.4]+2,'k')
text(12,2.25+2,'0.3 \muV')
text(90,4.5,'(C)','fontsize',20)

hold off

axis off
box off

ell = cab.ell;
xhat = 0.005*[1:4];
yhat = 0.01*[1:9];

I0 = stim.max*stim.tau*exp(1);
J0 = I0/(4*pi*sig)/lambda/sinh(ell/lambda);
y0 = stim.loc;

trueX = [J0 lambda y0];

save trapcabextraphi ell xhat yhat tim phi trueX

% finally, illustrate the important moments, V0 and V1

figure(4)
plot(x,V0/tau,'k','linewidth',1.5)
hold on
plot(x,V1/tau/tau,'r','linewidth',1.5)
legend('V_0/\tau','V_1/\tau^2','location','SE')
%text(stim.loc,max(V0)+7,'V_0/\tau','fontsize',14)
%text(stim.loc,max(V1)+7,'V_1/\tau^2','fontsize',14,'color','r')
text(.01,16,'(D)','fontsize',20)
xlabel('y  (cm)','fontsize',14)
ylabel('(mV)','fontsize',14)
hold off
box off

% check against analytical
% z = lambda;
% L = ell/z;  L0 = y0/z;
% tau = Cm/g_L;
% c = z*R_2/2/pi/cab.rad^2/sinh(L);
% I0 = stim.max*stim.tau*exp(1);
% I1 = 2*stim.max*(stim.tau)^2*exp(1);
% c1 = c*((2*I1+I0*tau)*cosh(L-L0)+(I0*tau/z)*(ell*cosh(L0)/sinh(L)+y0*sinh(L-L0)))
% c2 = c*cosh(L0)*(2*I1+I0*tau+(I0*tau/z)*(ell*cosh(L)/sinh(L)-y0*tanh(L0)));
% yL = 0:.001:y0;
% V1L = c1*cosh(yL/z) - I0*tau*R_2*cosh(L-L0)*yL.*sinh(yL/z)/2/pi/cab.rad^2/sinh(L);
% yR = y0:.001:ell;
% V1R = c2*cosh((ell-yR)/z) - I0*tau*R_2*cosh(L0)*(ell-yR).*sinh((ell-yR)/z)/2/pi/cab.rad^2/sinh(L);
% plot(yL,V1L)
% plot(yR,V1R)   % excellent fit so long as Tfin is long enough, say 40
% hold off

return

%  Recover cable and stimulus parameters from knowledge of extracelluar
%  potentials from trapcabextraphi.mat
%  generated by running trapcabextra.m
%
%   usage:  cabextrainv
%
function cabextrainv

close all
load trapcabextraphi

for j=1:9,
   Phi0(j,1) = trapz(tim',squeeze(phi(1,j,:)));       % zeroth moment of data
   Phi1(j,1) = trapz(tim',tim'.*squeeze(phi(1,j,:))); % first moment of data
end
    
[~,loc] = max(abs(Phi0));
Phi0rat = Phi0/Phi0(loc);

for j=1:6
    for k=1:99,
        
        y0 = k*0.001;
        lam = 0 + j*0.015;
        val = PhiMom([y0 lam],xhat(1),yhat,ell,Phi0rat,loc);
        E(k) = norm(val)^2;
        
    end
    
    plot(0.001*[1:99],E,'k')
    hold on
    
end

ylim([0 100])
xlabel('y_0  (cm)','fontsize',14)
ylabel('E(y_0,0.015m), m=1,...,6','fontsize',14)
text(0.09,90,'(A)','fontsize',20)
hold off

% now fix y0 and graph E wrt lam

for k=1:99,
    
    val = PhiMom([0.06 k*0.001],xhat(1),yhat,ell,Phi0rat,loc);
    E(k) = norm(val)^2;
    
end

figure(2)
plot(0.001*[1:99],E,'k')

ylim([0 0.1])
xlabel('\lambda  (cm)','fontsize',14)
ylabel('E(0.06,\lambda)','fontsize',14)
text(0.09,0.09,'(B)','fontsize',20)

% now solve for x = [tau; I1], 

lam = trueX(2);
J0 = trueX(1);
y0 = trueX(3);
X = trueX;

L = ell/lam;
L0 = y0/lam;
a1 = cosh(L-L0); a2 = cosh(L0);
b1 = (J0/2)*(L0*sinh(L-L0) + L*cosh(L0)/sinh(L) - cosh(L-L0));
b2 = (J0/2)*(L*cosh(L)*cosh(L0)/sinh(L) - cosh(L0) - L0*sinh(L0));
g1 = -(J0/2)*cosh(L-L0)/lam;
g2 = -(J0/2)*cosh(L0)/lam;

for i=1:9
    S1 = integral(@(y) ig1(y,X,xhat(1),yhat(i)),0,y0);
    S2 = integral(@(y) ig2(y,X,ell,xhat(1),yhat(i)),y0,ell);
    S3 = integral(@(y) ig3(y,X,xhat(1),yhat(i)),0,y0);
    S4 = integral(@(y) ig4(y,X,ell,xhat(1),yhat(i)),y0,ell);
    A(i,1) = a1*S1 + a2*S2 - lam*sinh(L)/sqrt(xhat(1)^2+(y0-yhat(i))^2);
    A(i,2) = b1*S1 + g1*S3 + b2*S2 + g2*S4;
end

x = A\Phi1;

J1 = x(1);
I1 = 4*pi*0.003*lam*sinh(L)*J1;
tau = x(2);

mc = A(:,2)./A(:,1);
rhs = Phi1./A(:,1);
figure(3)
plot(mc,rhs,'kx')   % raw data
hold on
LSA = ones(9,2);   % find best line fit through raw data
LSA(:,2) = mc;
x = LSA\rhs;
eta0 = -0.015;
eta1 = 0.02;
plot([eta0 eta1],x(1)+x(2)*[eta0 eta1],'k')
hold off
xlabel('\eta','fontsize',14)
ylabel('\zeta','fontsize',14)

return

function val = PhiMom(X,xhat,yhat,ell,Phi0rat,loc)
val = zeros(9,1);
for i=1:9,
    int1 = integral(@(y) ig1(y,X,xhat,yhat(i)),0,X(1));
    int2 = integral(@(y) ig2(y,X,ell,xhat,yhat(i)),X(1),ell);
    val(i) = (cosh((ell-X(1))/X(2))*int1 + cosh(X(1)/X(2))*int2 - ...
           X(2)*sinh(ell/X(2))/sqrt(xhat^2+(X(1)-yhat(i))^2));
end
val = val/val(loc);
val = val - Phi0rat;
return

function val = ig1(y,X,xhat,yhat)
val = cosh(y/X(2))./sqrt(xhat^2+(yhat-y).^2);
return

function val = ig2(y,X,ell,xhat,yhat)
val = cosh((ell-y)/X(2))./sqrt(xhat^2+(yhat-y).^2);
return

function val = ig3(y,X,xhat,yhat)
val = y.*sinh(y/X(2))./sqrt(xhat^2+(yhat-y).^2);
return

function val = ig4(y,X,ell,xhat,yhat)
val = (ell-y).*sinh((ell-y)/X(2))./sqrt(xhat^2+(yhat-y).^2);
return
