% Trapezoid on the fork with distributed input
% for a passive dendritive tree
%  usage:   trapforkd(dt,stim,pinc)
%
%               stim.t1 = start of current pulse (ms)
%               stim.t2 = end of current pulse (ms)
%               stim.amp = amplitude of current pulse (micro amps)
%               stim.Tfin = stopping time (ms)
% example:   stim = struct('t1',1,'t2',2,'amp',1e-6,'Tfin',8)
%               pinc = 4 
%

function trapforkd(dt,stim,pinc)

a = 1e-4*[1 1 1];
ell = [2.5 2.5 2.5]/100;
dx = .0001;
N = ell/dx;
A3 = 2*pi*a(3)*dx;
As = 4*pi*1e-6;
rho = A3/As;
R2 = 0.3;
gL = 1/15;
Cm = 1;
tau = Cm/gL;
lam = a/(2*R2*gL)/dx^2;   % lambda^2
r = a/a(3);
Hd = [2*lam(1)*ones(1,N(1)) 2*lam(2)*ones(1,N(2)) 2*lam(3)*ones(1,N(3)+1)];
Hd(1) = lam(1);
Hd(N(1)+1) = lam(2);
Hd(N(1)+N(2)+1) =  lam*r';
Hd(end) = rho*lam(3);
Hlen = length(Hd);

Hu = [-lam(1)*ones(1,N(1)-1) 0 -lam(2)*ones(1,N(2)) -lam(3)*ones(1,N(3))];
Hl = [-lam(1)*ones(1,N(1)-1) 0 -lam(2)*ones(1,N(2)-1) -r(2)*lam(2) -lam(3)*ones(1,N(3))];
Hl(end) = rho*Hl(end);

H = spdiags( [[Hl 0]' Hd' [0 Hu]'], -1:1, Hlen, Hlen);

H(N(1)+N(2)+1,N(1)) = -r(1)*lam(1);
H(N(1),N(1)+N(2)+1) = -lam(1);

I = speye(Hlen);

Bb = I+(I+H)*(dt/tau/2); 
[L,U] = lu(Bb); 

x3 = 0:dx:ell(3);
x1 = ell(3):dx:ell(3)+ell(1);
x2 = ell(3):dx:ell(3)+ell(2);
cx = cos(pi*x3'/2/ell(1));
xfoot = [cx(1:end-1); -cx(1:end-1); 0*cx];

v = zeros(Hlen,1);         % initial conditions
t = 0;

figure(2)
v1 = fliplr(v(1:N(1))');
v2 = fliplr(v(N(1)+1:N(1)+N(2))');
v3 = fliplr(v(N(1)+N(2)+1:end)');
plot3(x3,t*ones(size(x3)),v3,'r')
hold on
plot3(x2,t*ones(size(x2)),[v3(end) v2],'r')
plot3(x1,t*ones(size(x1)),[v3(end) v1],'k')

Nt = ceil(stim.Tfin/dt);

stim.amp = stim.amp*(dt/2)/A3/Cm;

f0 = stim.amp.*(t>stim.t1).*(t<stim.t2)*xfoot;
t = dt;
f1 = stim.amp.*(t>stim.t1).*(t<stim.t2)*xfoot;
figure(1)
f11 = fliplr(f1(1:N(1))');
f12 = fliplr(f1(N(1)+1:N(1)+N(2))');
f13 = fliplr(f1(N(1)+N(2)+1:end)');
plot3(x3,t*ones(size(x3)),f13,'r')
hold on
plot3(x2,t*ones(size(x2)),[f13(end) f12],'r')
plot3(x1,t*ones(size(x1)),[f13(end) f11],'k')

r = f0 + f1;

for j=2:Nt,

    v = U \ ( L \ r );

    if mod(j,pinc) == 0
         figure(2)
         v1 = fliplr(v(1:N(1))');
         v2 = fliplr(v(N(1)+1:N(1)+N(2))');
         v3 = fliplr(v(N(1)+N(2)+1:end)');
         plot3(x3,t*ones(size(x3)),v3,'r')
         plot3(x2,t*ones(size(x2)),[v3(end) v2],'r')
         plot3(x1,t*ones(size(x1)),[v3(end) v1],'k')

         figure(1)
         f11 = fliplr(f1(1:N(1))');
         f12 = fliplr(f1(N(1)+1:N(1)+N(2))');
         f13 = fliplr(f1(N(1)+N(2)+1:end)');
         plot3(x3,t*ones(size(x3)),f13,'r')
         hold on
         plot3(x2,t*ones(size(x2)),[f13(end) f12],'r')
         plot3(x1,t*ones(size(x1)),[f13(end) f11],'k')
    end

    t = j*dt;
    f0 = f1;
    f1 = stim.amp.*(t>stim.t1).*(t<stim.t2)*xfoot;

    r = 2*v - r + f0 + f1;

end

figure(1)
xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('I (\mu A)','fontsize',14)
axis tight
hold off

figure(2)
xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)
axis tight
hold off


