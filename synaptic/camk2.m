%  The dynamic response of the CaMKII system to a transient Ca2+ stimulus
%  w = camk2(dt,Tfin,castim)
%
%   dt = 0.025; 
%  Tfin = 6;  castim = struct('amp',100,'t1',1,'t2',.5)
%  Tfin = 60;  castim = struct('amp',10,'t1',10,'t2',5)
%

function w = camk2(dt,Tfin,castim)

close all

Nt = round(Tfin/dt);

pcam = struct('KM',.4,'KH1',4,'KH2',0.7,...
              'k1',5e-1,'k2',2,'k3',1,'k4',1e-3,...
              'P0',20,'D0',0.05,'I0',.1,...
              'nuCaN',1,'nuPKA',1);

w = [1 16/9 7/3 284/105 25/9 284/105 7/3 16/9 1]';

e = [0:10]';
B = spdiags([-e e],0:1,11,11);

P = zeros(11,1);

P(1) = pcam.P0;

A1 = (2/dt)*speye(11);
Fo = 0*A1;

Car = 0.1;
Ca3 = (Car/pcam.KH2)^3;
phi(2,1) = pcam.nuPKA*pcam.I0*(1+Ca3)/Ca3/pcam.nuCaN;
phi(1,1) = pcam.k4*pcam.D0/(pcam.k3*phi(2) + pcam.k4);
Phi1 = [-pcam.k4 0;-pcam.k4 -pcam.nuCaN*Ca3/(1+Ca3)];
V = [pcam.k4*pcam.D0  pcam.k4*pcam.D0+pcam.nuPKA*pcam.I0]';
U = -pcam.k3*phi(1)*phi(2)*[1 1]';

cat = zeros(Nt,1);
pat = cat;

for n=1:Nt,
    
    t = n*dt;
    
    Ca = castim.amp*(exp(-t/castim.t1)-exp(-t/castim.t2));
    
    Ca3 = (Ca/pcam.KH2)^3;
    Phi2 = [-pcam.k4 0;-pcam.k4 -pcam.nuCaN*Ca3/(1+Ca3)];
    
    phi = ((2/dt)*eye(2)-Phi2) \ ( ((2/dt)*eye(2)+Phi1)*phi + 2*(U+V) );
    
    U = -pcam.k3*phi(1)*phi(2)*[1 1]';
    Phi1 = Phi2;
    
    b1 = pcam.k2*phi(1)/(pcam.KM + sum([1:10]'.*P(2:11)));
    
    ct = (Ca/pcam.KH1)^4;
    hct = ct/(1+ct);
    f = pcam.k1*hct*ones(10,1);

   f(1) = 10*hct*f(1);
   f(2:10) = w(1:9).*f(2:10);
   f(11) = 0;
   Fn = spdiags([f -f],-1:0,11,11);
    
    P = (A1-Fn-b1*B)\((A1+Fo+b1*B)*P);

    cat(n) = Ca;
    pat(n) = sum([1:10]'.*P(2:11));
    
    Fo = Fn;
    
end

t = linspace(0,Tfin,Nt)';
subplot(2,1,1)
plot(t,cat,'k')
box off
ylabel('[Ca^{2+}]  (\muM)','fontsize',14)
subplot(2,1,2)
plot(t,pat,'k')
box off
xlabel('t  (s)','fontsize',14)
ylabel('[P]_A  (\muM)','fontsize',14)

