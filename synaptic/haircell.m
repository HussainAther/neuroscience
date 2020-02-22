% The hair cell of Hudspeth and Lewis

function haircell1

haircell(.01,90,20,'r')
hold on
haircell(.01,90,80,'k')
hold off
legend('I_{stim}=20 pA','I_{stim}=80 pA')

return

function haircell(dt,tfin,I0,col)

Nt = tfin/dt + 1;
V = zeros(1,Nt);t = V;Istim = V;
Cm = 15;  % pF
g = struct('Ca', 4.14,'KCa', 16.8,'L', 1);   % nS
E = struct('Ca',100,'K',-80,'L', -30);  
D = [0.2 0 0.2];
K0 = [6 45 20]*1e-3;    %  mM
kbac = [0.3 5 1.5];   % 1/ms
F = 9.6485e4; % coul/mol
gamma = 0.02/(2*F*1.25*3.4e-5);
p = struct('gamma',gamma,'KS', 2.800, 'betaC', 1, 'alphaC0', .45, 'Va', 33, 'VT', 25.8);

Vr = fsolve(@(v) irest(v,g,E,D,K0,p),-70);
v = Vr;
m = am(v)/(am(v)+bm(v));
ICa = g.Ca*m^3*(v-E.Ca);
c = -p.gamma*ICa/p.KS;
rat = c.*exp(-D*2*v/p.VT)./K0;
R = p.betaC*exp(v/p.Va)/p.alphaC0;
C0 = 1./(1+rat(1)*(1+rat(2)+(rat(3)+1)*rat(2)*R));
C1 = rat(1)*C0;
C2 = rat(2)*C1;
O2 = R*C2;
O3 = rat(3)*O2;
s = [C0 C1 C2 O2 O3]';

if 1==0
   v = -80:.1:0;
   [ICa, IKCa, IL] =  icurr(v,g,E,D,K0,p);
   plot(v,ICa,v,IKCa,v,IL)
   legend('ICa','IKCa','IL')
   hold on
   plot(v,ICa+IKCa+IL,'k')
   hold off
end

I5 = eye(5);
V(1) = v;

A = zeros(5);
A(1,2) = kbac(1);
A(2,3) = kbac(2);
A(3,3) = -(kbac(2)+p.betaC);
A(4,[3  5]) = [p.betaC kbac(3)];
A(5,5) = -kbac(3);

for j=2:Nt

    a = am(v); apb = a + bm(v);
    m = (m + dt*a)/(1 + dt*apb);
    ICa = g.Ca*m^3*(v-E.Ca);

    c = (c - dt*p.gamma*ICa)/(1+dt*p.KS);

    K = K0.*exp(D*2*v/p.VT);
    k = kbac./K;
    alphaC = p.alphaC0*exp(-v/p.Va);

    A(1,1) = -k(1)*c;
    A(2,1) = k(1)*c;
    A(2,2) = -(kbac(1)+k(2)*c);
    A(3,2) = k(2)*c;
    A(3,4) = alphaC;
    A(4,4) = -(alphaC + k(3)*c);
    A(5,4) = k(3)*c;

    s = (I5 - dt*A)\s;

    O23 = s(4) + s(5);

    IKCa = g.KCa*O23.*(v-E.K);

    t(j) = dt*(j-1);
    
    Istim(j) = (t(j)>5)*(t(j)<55)*I0;
    
    top = v*Cm/dt + g.Ca*m^3*E.Ca + g.KCa*O23*E.K + g.L*E.L + Istim(j);
    bot = Cm/dt + g.Ca*m^3 + g.KCa*O23 + g.L;

    v = top/bot;
    V(j) = v;
    
end

plot(t,V,col)
ylabel('V  (mV)','fontsize',14)
xlabel('t  (ms)','fontsize',14)
set(gca,'tickdir','out')
box off

return

function [ICa, IKCa, IL] = icurr(v,g,E,D,K0,p)
m = am(v)./(am(v)+bm(v));
ICa = g.Ca*m.^3.*(v-E.Ca);
c = -p.gamma*ICa/p.KS;
for j=1:3
  rat(j,:) = c.*exp(-D(j)*2*v/p.VT)/K0(j);
end
R = p.betaC*exp(v/p.Va)/p.alphaC0;
C0 = 1./(1+rat(1,:).*(1+rat(2,:)+(rat(3,:)+1).*rat(2,:).*R));
O23 = C0.*(rat(3,:)+1).*rat(2,:).*rat(1,:).*R;
IKCa = g.KCa*O23.*(v-E.K);
IL = g.L*(v-E.L);
val = ICa + IKCa + IL;

function val = irest(v,g,E,D,K0,p)
m = am(v)/(am(v)+bm(v));
ICa = g.Ca*m^3*(v-E.Ca);
c = -p.gamma*ICa/p.KS;
rat = c.*exp(-D*2*v/p.VT)./K0;
R = p.betaC*exp(v/p.Va)/p.alphaC0;
C0 = 1./(1+rat(1)*(1+rat(2)+(rat(3)+1)*rat(2)*R));
O23 = C0*(rat(3)+1)*rat(2)*rat(1)*R;
IKCa = g.KCa*O23*(v-E.K);
IL = g.L*(v-E.L);
val = ICa + IKCa + IL;

% ICa functionals

function val = am(v)
val = 0.00097*exp((v+70)/6.17) + .940;

function val = bm(v)
val = 22.8*exp(-(v+70)/8.01) + .510;

% Drive the hair cell of Hudspeth and Lewis with hair displacement
% over a range of frequencies
function haircell2

for j=1:100

    f(j) = j*.002;
    v = haircell(0.005,800,f(j));
    PP(j) = max(v(20000:end))-min(v(20000:end));

end

plot(f*1000,PP,'k')
xlabel('f  (Hz)','fontsize',16)
ylabel('P  (mV)','fontsize',16)
set(gca,'tickdir','out')
box off

return

function V = haircell(dt,tfin,f)

Nt = tfin/dt + 1;
V = zeros(1,Nt);t = V;Istim = V;INa = V;IK = V;IA = V;
Cm = 15;  % pF
g = struct('Ca', 4.14,'KCa', 16.8,'L', 1);   % nS
E = struct('Ca',100,'K',-80,'L', -30);  
D = [0.2 0 0.2];
K0 = [6 45 20]*1e-3;    %  mM
kbac = [0.3 5 1.5];   % 1/ms
F = 9.6485e4; % coul/mol
gamma = 0.02/(2*F*1.25*3.4e-5);
%pause
p = struct('gamma',gamma,'KS', 2.800, 'betaC', 1, 'alphaC0', .45, 'Va', 33, 'VT', 25.8);

Vr = fsolve(@(v) irest(v,g,E,D,K0,p),-70);
v = Vr;
m = am(v)/(am(v)+bm(v));
ICa = g.Ca*m^3*(v-E.Ca);
c = -p.gamma*ICa/p.KS;
rat = c.*exp(-D*2*v/p.VT)./K0;
R = p.betaC*exp(v/p.Va)/p.alphaC0;
C0 = 1./(1+rat(1)*(1+rat(2)+(rat(3)+1)*rat(2)*R));
C1 = rat(1)*C0;
C2 = rat(2)*C1;
O2 = R*C2;
O3 = rat(3)*O2;
s = [C0 C1 C2 O2 O3]';

if 1==0
   v = -80:.1:0;
   [ICa, IKCa, IL] =  icurr(v,g,E,D,K0,p);
   plot(v,ICa,v,IKCa,v,IL)
   legend('ICa','IKCa','IL')
   hold on
   plot(v,ICa+IKCa+IL,'k')
   hold off
end

I5 = eye(5);
V(1) = v;

A = zeros(5);
A(1,2) = kbac(1);
A(2,3) = kbac(2);
A(3,3) = -(kbac(2)+p.betaC);
A(4,[3  5]) = [p.betaC kbac(3)];
A(5,5) = -kbac(3);

for j=2:Nt

    a = am(v); apb = a + bm(v);
    m = (m + dt*a)/(1 + dt*apb);
    ICa = g.Ca*m^3*(v-E.Ca);

    c = (c - dt*p.gamma*ICa)/(1+dt*p.KS);

    K = K0.*exp(D*2*v/p.VT);
    k = kbac./K;
    alphaC = p.alphaC0*exp(-v/p.Va);

    A(1,1) = -k(1)*c;
    A(2,1) = k(1)*c;
    A(2,2) = -(kbac(1)+k(2)*c);
    A(3,2) = k(2)*c;
    A(3,4) = alphaC;
    A(4,4) = -(alphaC + k(3)*c);
    A(5,4) = k(3)*c;

    s = (I5 - dt*A)\s;

    O23 = s(4) + s(5);

    IKCa = g.KCa*O23.*(v-E.K);

    t(j) = dt*(j-1);
    
    x = 20*sin(f*2*pi*t(j));     % nm
    RT = 8.314*298;
    eA = exp((3140 - 41.9*x)/RT);
    eB = exp((1050 - 8.4*x)/RT);
    gT = 3/(1+eB*(1+eA));
    
    top = v*Cm/dt + g.Ca*m^3*E.Ca + g.KCa*O23*E.K + g.L*E.L - gT*v;
    bot = Cm/dt + g.Ca*m^3 + g.KCa*O23 + g.L;

    v = top/bot;
    V(j) = v;
    
end

plot(t,V)
ylabel('V  (mV)','fontsize',14)
xlabel('t  (ms)','fontsize',14)
set(gca,'tickdir','out')
box off
drawnow

return

function [ICa, IKCa, IL] = icurr(v,g,E,D,K0,p)
m = am(v)./(am(v)+bm(v));
ICa = g.Ca*m.^3.*(v-E.Ca);
c = -p.gamma*ICa/p.KS;
for j=1:3
  rat(j,:) = c.*exp(-D(j)*2*v/p.VT)/K0(j);
end
R = p.betaC*exp(v/p.Va)/p.alphaC0;
C0 = 1./(1+rat(1,:).*(1+rat(2,:)+(rat(3,:)+1).*rat(2,:).*R));
O23 = C0.*(rat(3,:)+1).*rat(2,:).*rat(1,:).*R;
IKCa = g.KCa*O23.*(v-E.K);
IL = g.L*(v-E.L);
val = ICa + IKCa + IL;

function val = irest(v,g,E,D,K0,p)
m = am(v)/(am(v)+bm(v));
ICa = g.Ca*m^3*(v-E.Ca);
c = -p.gamma*ICa/p.KS;
rat = c.*exp(-D*2*v/p.VT)./K0;
R = p.betaC*exp(v/p.Va)/p.alphaC0;
C0 = 1./(1+rat(1)*(1+rat(2)+(rat(3)+1)*rat(2)*R));
O23 = C0*(rat(3)+1)*rat(2)*rat(1)*R;
IKCa = g.KCa*O23*(v-E.K);
IL = g.L*(v-E.L);
val = ICa + IKCa + IL;

% ICa functionals

function val = am(v)
val = 0.00097*exp((v+70)/6.17) + .940;

function val = bm(v)
val = 22.8*exp(-(v+70)/8.01) + .510;
