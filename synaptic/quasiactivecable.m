% Solve the quasi-active cable problem via staggered Euler 
% to accomodate observation of the potential at both x = l/10 and
% x = 9l/10
%
% usage:	stEQcabBT2(cab,stim,pinc)
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
% examples:	
%
%  cab = struct('rad',1e-4,'ell',.1,'dx',1e-3,'dt',0.05)
%  stim = struct('t1',1,'t2',2,'amp',5e-6,'loc',0.2,'Tfin',8)
%  pinc = 10 
%
%  N = 50; t1rand = rand(1,N)*N; locrand = rand(N,1)*cab.ell;
%  stim = struct('t1',t1rand,'t2',t1rand+1,'amp',5e-5,'loc',locrand,'Tfin',max(t1rand)+4)
%

function stEQcabBT2(cab,stim,pinc)

Cm = 1;		% micro F / cm^2
R2 = 0.3;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;		% patch length
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints

E = struct('K', -77, 'Na', 56, 'Cl', -68);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
%g = struct('K', 0, 'Na', 0, 'Cl', 0.3);

eloc = max(round(Nx*stim.loc/cab.ell),1);
obloc = [round(Nx*0.1) round(Nx*0.9)]          % observer

e = ones(Nx,1);
S = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = 1/dx/dx;
S(Nx,Nx) = 1/dx/dx;
S = (cab.rad/2/R2)*S;
Z = 0*speye(Nx);

options = optimset('Jacobian','on');
V = fsolve(@(V) Iss(V,E,g,S),-70*e,options);    % initial conditions

a = am(V); b = bm(V);
taum = 1./(a+b);
m = a.*taum;
dm = (dam(V).*b - a.*dbm(V)).*taum.^2;

a = ah(V); b = bh(V);
tauh = 1./(a+b);
h = a.*tauh;
dh = (dah(V).*b - a.*dbh(V)).*tauh.^2;

a = an(V); b = bn(V);
taun = 1./(a+b);
n = a.*taun;
dn = (dan(V).*b - a.*dbn(V)).*taun.^2;

Q11 = -(S + spdiags(g.Cl+g.Na*m.^3.*h+g.K*n.^4,0,Nx,Nx))/Cm;
Q12 = spdiags(3*g.Na*m.^2.*h*E.Na/Cm,0,Nx,Nx);
Q13 = spdiags(g.Na*m.^3*E.Na/Cm,0,Nx,Nx);
Q14 = spdiags(4*g.K*n.^3.*E.K/Cm,0,Nx,Nx);

Q21 = spdiags(dm./taum,0,Nx,Nx);
Q22 = spdiags(-1./taum,0,Nx,Nx);

Q31 = spdiags(dh./tauh,0,Nx,Nx);
Q33 = spdiags(-1./tauh,0,Nx,Nx);

Q41 = spdiags(dn./taun,0,Nx,Nx);
Q44 = spdiags(-1./taun,0,Nx,Nx);

B = [Q11 Q12 Q13 Q14;
     Q21 Q22  Z   Z;
     Q31  Z  Q33  Z;
     Q41  Z   Z  Q44];

[L,U] = lu(speye(4*Nx)-B*dt/2);

Nt = ceil(stim.Tfin/dt);

y = zeros(4*Nx,1);
yob = zeros(2,Nt);

t = 0;
figure(1)
plot3(x,t*e,y(1:Nx),'k')
hold on

f0 = stim.amp*(t>stim.t1).*(t<stim.t2)/A;
t = dt;
f1 = stim.amp*(t>stim.t1).*(t<stim.t2)/A;

rhs = (speye(4*Nx)+B*dt/2)*y; 
rhs(eloc) = rhs(eloc) + (f0+f1)'*dt/2;

for j=2:Nt,
    
    y = U \ ( L \ rhs );
    
    yob(:,j) = y(obloc);

    if mod(j,pinc) == 0
        plot3(x,t*e,y(1:Nx),'k')
    end

    t = j*dt;
    f0 = f1;
    f1 = stim.amp*(t>stim.t1).*(t<stim.t2)/A;

    rhs = 2*y - rhs;
    rhs(eloc) = rhs(eloc) + (f0+f1)'*dt/2;

end

xlabel('x  (cm)','fontsize',14)
ylabel('t  (ms)','fontsize',14)
zlabel('v  (mV)','fontsize',14)

hold off

figure(2)
tim = linspace(0,stim.Tfin,Nt);
%plot(tim,yob,'k')
plot(tim,yob(1,:),'k')
hold on
plot(tim,yob(2,:),'r')

% Construct the reduced system

e = ones(4*Nx,1); e(Nx+1:4*Nx) = 0;
C = spdiags(e,0,4*Nx,4*Nx);

%D = zeros(1,4*Nx); D(obloc) = 1;
D = zeros(2,4*Nx); D(1,obloc(1)) = 1; D(2,obloc(2)) = 1;

B = full(B);

RP = lyapchol(B,full(C));
U = RP';

RQ = lyapchol(B',D');
L = RQ';

[Y,Sig,X] = svd(U'*L);

sig = diag(Sig);

figure(3)
semilogy(sig(1:100),'kx-')
box off
xlabel('index','fontsize',14)
ylabel('singular value','fontsize',14)

sigmh = diag(1./sqrt(sig));
Phi = sigmh*X'*L';
PhiI = U*Y*sigmh;

k = 5;
Bhat = Phi*B*PhiI;
BhatRed = Bhat(1:k,1:k);
Chat = Phi*C;
ChatRed = Chat(1:k,:);
Dhat = D*PhiI;
DhatRed = Dhat(:,1:k);

% solve the reduced system

[L,U] = lu(eye(k)-BhatRed*dt/2);

xi = zeros(k,1);
yobRed = zeros(2,Nt);
u = zeros(4*Nx,1);

t = 0;
f0 = stim.amp*(t>stim.t1).*(t<stim.t2)/A;
t = dt;
f1 = stim.amp*(t>stim.t1).*(t<stim.t2)/A;
u(eloc) = (f0+f1)'*dt/2;
rhs = (eye(k)+BhatRed*dt/2)*xi + ChatRed*u;

for j=2:Nt,
    
    xi = U \ ( L \ rhs );
    
    yobRed(:,j) = DhatRed*xi;

    t = j*dt;
    f0 = f1;
    f1 = stim.amp*(t>stim.t1).*(t<stim.t2)/A;
    u(eloc) = (f0+f1)'*dt/2;
    rhs = 2*xi - rhs + ChatRed*u;

end

figure(2)
plot(tim,yobRed(1,:),'k--')
plot(tim,yobRed(2,:),'r--')
hold off
box off
axis tight
xlabel('t  (ms)','fontsize',14)
ylabel('(mV)','fontsize',14)
%legend('full, N=400','reduced, k=5','location','best')

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

