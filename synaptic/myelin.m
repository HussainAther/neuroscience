% Solve the myelinated (myelin) active cable problem via staggered Euler
% and return the wavespeed
%
% usage:  spd = myelins(cab,stim,pinc)
%
%       cab.rad = fiber radius (cm)
%	cab.id = internodal distance (microns)
%	cab.nnor = number of nodes of Ranvier
%       cab.dt = timestep (ms)
%
%	stim.t1 = start of current pulse (ms)
%	stim.t2 = end of current pulse (ms)
%       stim.amp = amplitude of current pulse (micro amps)
%       stim.loc = location of current pulse (cm)
%       stim.Tfin = stopping time (ms)
%
%       pinc = number of time steps between plots
%
%       computed wave speed
% e.g.:	
%       cab = struct('rad',2e-4,'id',2000,'nnor',10,'dt',.05)
%       stim = struct('t1',1,'t2',2,'amp',400e-6,'loc',1,'Tfin',15)
%       pinc = 10
%
%       spd = myelins(cab,stim,pinc)
%

function spd = myelins(cab,stim,pinc)

dx = 1e-4;				% space step
Nx = 2*cab.nnor + (cab.nnor-1)*cab.id;

pat = zeros(Nx,1);
pat(1:cab.id+2:Nx) = 1;
pat(2:cab.id+2:Nx) = 1;

R2 = 0.3;             % k Ohm cm
dt = cab.dt;
A = 2*pi*cab.rad*dx;		% patch surface area
ell = Nx*dx;
x = dx/2:dx:ell-dx/2;	% vector of patch midpoints

Nt = ceil(stim.Tfin/dt);  % number of time steps

E = struct('K', -77, 'Na', 56, 'Cl', -68);

e = ones(Nx,1);

g.K = 36*pat;
g.Na = 120*pat;         	

Cm = 0.01*e + (1-0.01)*pat;
g.Cl = 0.003*e + (0.3-0.003)*pat;

Iapp = stim.amp/A;
eloc = stim.loc;

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
Va(1) = V(cab.id+2);
Vb(1) = V(Nx-cab.id-1);

t = 0;
if pinc 
   figure(1)
   plot(x,V,'k')
   xlim([0 ell])
   box off
   xlabel('x  (cm)')
   ylabel('V_{rest}  (mV)')
   figure(2)
   plot3(x,t*ones(Nx,1),V,'color','k')
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

      d = g.Na.*m3h + g.K.*n4 + g.Cl;

      f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl.*E.Cl;
      f(eloc) = f(eloc) + I;

      B(1:Nx+1:end) = dB + d + 2*Cm./dt;         % update the diagonal

      Vmid = B\(2*Cm.*V/dt + f);
      
      V = 2*Vmid - V;
      
      Va(j) = V(cab.id+2);
      Vb(j) = V(Nx-cab.id-1);

      if mod(j,pinc) == 0
         plot3(x,t*ones(Nx,1),V,'color','k')
      end

end

[val,inda] = max(Va);
[val,indb] = max(Vb);
     
spd = (x(Nx-cab.id-1)-x(cab.id+2))/(((indb-1)*cab.dt)-((inda-1)*cab.dt));

if pinc
    xlabel('x  (cm)','fontsize',14)
    ylabel('t  (ms)','fontsize',14)
    zlabel('V  (mV)','fontsize',14)
    box off
    hold off
end

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

val = B*V + g.Na.*m3h.*(V-E.Na) + g.K.*n4.*(V-E.K) + g.Cl.*(V-E.Cl);

dj = g.Na.*((3*dm.*m.^2.*h + m.^3.*dh).*(V-E.Na) + m3h) + ...
     g.K.*(4*dn.*n.^3.*(V-E.K) + n4) + g.Cl;
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

