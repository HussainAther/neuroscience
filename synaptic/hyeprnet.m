% Solve the Pinsky-Rinzel 2 compartment random (or prescribed) CA3 
% net with N cells connected with density den
%
%  usage   A = hyEprnet(dt,Tfin,N,den,gAMPA,gNMDA,TAMPA,adj)

function A = hyEprnet(dt,Tfin,N,den,gAMPA,gNMDA,TAMPA,adj)

Nt = ceil(Tfin/dt);

if nargin == 8
    A = adj;
else
    A = ceil(sprand(N,N,den));
    A(1:N+1:end) = 0;
end
figure(1)
hist(sum(A,2))

%N = 30;		Simple ring
%A = zeros(N);
%A(1,N) = 1;
%for i=2:N,
%    A(i,i-1) = 1;
%end

e = ones(N,1);
Vs = -4.6*e;  Vd = -4.5*e;
h = 0.999*e;  n = 0.001*e;
s = 0.009*e;  r = 0.007*e;
q = 0.01*e;   c = 0.2*e;
x = 0*e;      y = 0*e;

g = struct('L', 0.1,'Na', 30, 'K', 15, 'Ca', 10, 'AHP', 0.8, 'KC', 15,...
           'AMPA', gAMPA, 'NMDA', gNMDA);
E = struct('Na', 120, 'Ca', 140, 'K', -15, 'L', 0, 'syn', 60);

gc = 2.1;
p = 0.5;
Cm = 3;
oodt = 1/dt;

spk = zeros(1,Nt);
tic

hi = zeros(2*N,1); lo = hi; tdi = hi; rhs = hi;
hi(2:2:2*N) = -dt*gc/p;
lo(1:2:2*N) = -dt*gc/(1-p);
di = (Cm + dt*(g.L+gc/p))*ones(2*N,1);
A2 = spdiags([lo di hi], -1:1, 2*N, 2*N);

odx = 1/(1+dt/2);  ody = 1/(1+dt/150);
dtA = dt*A;

for j=2:Nt,

    a = ah(Vs);  abm = (a+bh(Vs))/2;
    h = ((oodt-abm).*h + a)./(oodt + abm); 

    a = an(Vs);  abm = (a+bn(Vs))/2;
    n = ((oodt-abm).*n + a)./(oodt+abm);

    a = am(Vs);  
    m = a./(a + bm(Vs));

    a = as(Vd);  abm = (a + bs(Vd))/2;
    s = ((oodt - abm).*s + a)./(oodt + abm);

    a = ar(Vd);  abm = (a + br(Vd))/2;
    r = ((oodt - abm).*r + a)./(oodt + abm);

    a = aq(c);  abm = (a + 0.001)/2;
    q = ((oodt - abm).*q + a)./(oodt + abm);

    chi = min(c/250,1);
    
    t = (j-1)*dt;

    e20 = Vs>20;
    x = (x + dtA*e20)*odx;
    spk(j) = sum(e20);

    e10 = Vs>10;
    y = (y + dtA*e10)*ody;
    y = min(y,125);

    gsyn = (g.AMPA*(t<TAMPA)*x + g.NMDA*y.*M(Vd))/(1-p);

    tNa = g.Na*m.^2.*h;
    tK = g.K*n;
    tCa = g.Ca*s.^2;
    tCo = g.KC*chi.*r + g.AHP*q;

    tdi(1:2:2*N) = dt*(tNa + tK);
    tdi(2:2:2*N) = dt*(tCa + tCo + gsyn);

    A2(1:2*N+1:end) = di + tdi;

    rhs(1:2:2*N) = Cm*Vs + dt*(tNa*E.Na + tK*E.K);
    rhs(1) = rhs(1) + dt*(-0.5+30*(t>1)*(t<3))/p;
    rhs(2:2:2*N) = Cm*Vd + dt*(tCa*E.Ca + tCo*E.K + gsyn*E.syn);

    V = A2 \ rhs;

    Vs = V(1:2:2*N);  Vd = V(2:2:2*N);

    ICa = g.Ca*(s.^2).*(Vd-E.Ca);

    c = (c - dt*0.13*ICa)/(1 + dt*0.075);

end

toc
t = linspace(0,Tfin,Nt);
figure(2)
plot(t,spk/N,'k')
ylim([-.1 1.1])
box off
set(gca,'tickdir','out')
xlabel('t  (ms)','fontsize',14)
ylabel('Burst Fraction','fontsize',14)

return

function val = am(v)
val = 0.32*(13.1-v)./(exp((13.1-v)/4)-1);

function val = bm(v)
val = 0.28*(v-40.1)./(exp((v-40.1)/5)-1);

function val = an(v)
val = 0.016*(35.1-v)./(exp((35.1-v)/5)-1);

function val = bn(v)
val = 0.25*exp(0.5-0.025*v);

function val = ah(v)
val = 0.128*exp((17-v)/18);

function val = bh(v)
val = 4./(1+exp((40-v)/5));

function val = as(v)
val = 1.6./(1+exp(-0.072*(v-65)));

function val = bs(v)
val = 0.02*(v-51.1)./(exp((v-51.1)/5)-1);

function val = ar(v)
val = (v<50).*exp((v-10)/11-(v-6.5)/27)/18.975 + (v>50)*2.*exp((6.5-v)/27);

function val = br(v)
val = (v<50).*(2*exp((6.5-v)/27) - exp((v-10)/11-(v-6.5)/27)/18.975);

function val = aq(Ca)
val = min((0.00002)*Ca,0.01);

function val = bq(v)
val = 0.001;

function val = M(v)
val = 1./(1+0.28*exp(-0.062*(v-60)));
