% Determine wave speed. Investigate the impact of cable radius and gNa on 
% the velocity of the action potential that propagates down the active uniform cable.

function stEcabsdriver

cab = struct('rad',1e-3,'ell',1,'dx',5e-3,'dt',.01);
stim = struct('t1',1,'t2',2,'amp',5e-2,'loc',5e-3,'Tfin',20);

for j=1:45
    arad(j) = 5e-5 + (j-1)*1e-5;
    cab.rad = arad(j);
    spd(j) = stEcabs2(cab,stim,0);
end

figure(1)
plot(arad*1e4,spd*10,'kx-')
ylabel('Wave Speed (m/s)','fontsize',14)
xlabel('Cable Radius (\mum)','fontsize',14)
box off
hold off

%
% stEcabs.m
%
% Cox & Gabbiani, Mathematics for Neuroscientists
%
% determine wave speed in stEcab
%
% usage:	stEcab(cab,stim,pinc)
%
%		cab.rad = fiber radius (cm)
%		cab.ell = fiber length (cm)
%       cab.dx = space step (cm)
%       cab.dt = timestep (ms)
%		stim.t1 = start of current pulse (ms)
%		stim.t2 = end of current pulse (ms)
%       stim.amp = amplitude of current pulse (micro amps)
%       stim.loc = location of current pulse (cm)
%       stim.Tfin = stopping time (ms)
%       pinc = number of time steps between plots
%
% example:	cab = struct('rad',1e-3,'ell',1,'dx',1e-2,'dt',0.1)
%           stim = struct('t1',1,'t2',2,'amp',5e-2,'loc',1e-2,'Tfin',8)
%

function spd = stEcabs2(cab,stim,pinc)

Cm = 1;		% micro F / cm^2
R2 = 0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;		% patch length
A = 2*pi*cab.rad.*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

E = struct('K', -77, 'Na', 56, 'Cl', -68);
g = struct('K', 36, 'Na', 120, 'Cl', 0.3);

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
e1(eloc) = 1;
Iapp = stim.amp*e1/A;

e = ones(Nx,1);
B = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
B(1,1) = 1/dx/dx;
B(Nx,Nx) = 1/dx/dx;
B = (cab.rad/2/R2)*B;
dB = diag(B);

options = optimset('Jacobian','on');
V = fsolve(@(V) Iss(V,E,g,B),-70*e,options);    % initial conditions

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

Va = zeros(Nt,1);
Vb = zeros(Nt,1);
Va(1) = V(10);
Vb(1) = V(Nx-10);

figure(1)

t= 0;
if pinc > 0
   figure(2)
   plot3(x,t*ones(Nx,1),V)
   hold on
end

for j=2:Nt,
           
      t = (j-1)*dt;

      I = Iapp*(t-dt/2>stim.t1)*(t-dt/2<stim.t2);

      a = an(V);  c = (a+bn(V))/2;
      n = ( (1/dt-c).*n + a) ./ (1/dt + c); n4 = n.^4;

      a = am(V);  c = (a+bm(V))/2;
      m = ( (1/dt-c).*m + a) ./ (1/dt + c);

      a = ah(V);  c = (a+bh(V))/2;
      h = ( (1/dt-c).*h + a) ./ (1/dt + c); m3h = m.^3.*h;

      d = g.Na*m3h + g.K*n4 + g.Cl;

      f = g.Na*m3h*E.Na + g.K*n4*E.K + g.Cl*E.Cl + I;

      B(1:Nx+1:end) = dB + d + 2*Cm/dt;         % update the diagonal

      Vmid = B\(2*Cm*V/dt + f);
      
      V = 2*Vmid - V;
      
      Va(j) = V(10);
      Vb(j) = V(Nx-10);
      
      if mod(j,pinc) == 0
         figure(2)
         plot3(x,t*ones(Nx,1),V,'color','k')
         xlabel('Distance (m)','fontsize',14)
         ylabel('Time (s)','fontsize',14)
         zlabel('Voltage(mv)','fontsize',14)
      end
      
end
      
     [val,inda] = max(Va);
     [val,indb] = max(Vb);
     
     spd = (x(Nx-10)-x(10))/(((indb-1)*cab.dt)-((inda-1)*cab.dt));

     hold off
     
return

function [val, jac] = Iss(V,E,g,B)
Nx = length(V);
a = am(V); b = bm(V);
m = a./(a+b);
dm = (dam(V).*b - a.*dbm(V))./(a+b).^2;

a = ah(V); b = bh(V);
h = a./(a+b);
dh = (dah(V).*b - a.*dbh(V))./(a+b).^2;

a = an(V); b = bn(V);
n = a./(a+b);
dn = (dan(V).*b - a.*dbn(V))./(a+b).^2;

m3h = m.^3.*h;
n4 = n.^4;

val = B*V + g.Na*m3h.*(V-E.Na) + g.K*n4.*(V-E.K) + g.Cl.*(V-E.Cl);

dj = g.Na*((3*dm.*m.^2.*h + m.^3.*dh).*(V-E.Na) + m3h) + ...
     g.K*(4*dn.*n.^3.*(V-E.K) + n4) + g.Cl;
jac = B + spdiags(dj,0,Nx,Nx);

function val = an(v)
val = .01*(10-(v+71))./(exp(1-(v+71)/10)-1);

function val = dan(v)
tmp = exp(-(61+v)/10);
val = -( tmp.*(71+v) - 10 )./(tmp-1).^2/1000;

function val = bn(v)
val = .125*exp(-(v+71)/80);

function val = dbn(v)
val = -exp(-(v+71)/80)/640;

function val = am(v)
val = .1*(25-(v+71))./(exp(2.5-(v+71)/10)-1);

function val = dam(v)
tmp = exp(-(46+v)/10);
val = -( tmp.*(56+v) - 10 )./(tmp-1).^2/100;

function val = bm(v)
val = 4*exp(-(v+71)/18);

function val = dbm(v)
val = -(2/9)*exp(-(v+71)/18);

function val = ah(v)
val = 0.07*exp(-(v+71)/20);

function val = dah(v)
val = -(7/2000)*exp(-(v+71)/20);

function val = bh(v)
val = 1./(exp(3-(v+71)/10)+1);

function val = dbh(v)
tmp = exp(-(v+41)/10);
val = tmp./(tmp+1).^2/10;

