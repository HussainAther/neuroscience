% Solve the active cable via sparse hybrid trapezoid
% with 3 calcium channels AND NaCa exchange
% and solution of the u=[c;b] reaction diffusion system
%
%  and gKCa
%
% usage: [vmid, t] = hyEcabCa3train(cab,stim,pinc)
%
%		cab.rad = cable radius (cm)
%		cab.ell = cable length (cm)
%       cab.N = number of compartments
%       cab.dt = timestep (ms)
%		stim.pw = start of current pulse (ms)
%		stim.per = end of current pulse (ms)
%       stim.Tfin = stopping time (ms)
%		stim.Iapp = size of current pulse (micro amps)
%       stim.BT = total buffer concentration
%       stim.knaca = NaCa rate
%       stim.Jpmca = pmca pump rate
%       pinc = number of time steps between plots
%
% example:	
%           cab = struct('rad',.0001,'ell',.1,'N',200,'dt',.02)
%           stim = struct('pw',1,'per',20,'Tfin',100,'Iapp',3e-4,...
%           'BT',5e2,'knaca',100,'Jpmca',1e-5,'gKCa',10)
%           pinc = 20
%           [vmid, t] = hyEcabCa3train(cab,stim,pinc)
%

function [vmid, t] = hyEcabCa3train(cab,stim,pinc)

Cm = 1;         % micro F / cm^2
R2 = 0.3;       % 0.034;  % k Ohm cm
N = cab.N;
ell = cab.ell;
dx = ell/N;		% patch length
A = 2*pi*cab.rad*dx;		% patch surface area
x = dx/2:dx:cab.ell-dx/2;	% vector of patch midpoints
Nmid = round(N/2);

Nt = ceil(stim.Tfin/cab.dt)+1;
ICaT = zeros(Nt+1,1); ICaN = ICaT; ICaL = ICaT; vmid = ICaT;
ICa = ICaT; ca = ICa; bca = ICa; INaCa = ca; IK = ICa;

E = struct('K', -77, 'Na', 56, 'Cl', -68);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15, ...
           'CaT', .25, 'CaN', 2.5, 'CaL', 2.5,'KCa',stim.gKCa);
        
e1 = zeros(N,1);
e1(Nmid/2) = 1;
Iapp = stim.Iapp*e1/A;
 
e = ones(N,1);

Na3 = (50/440)^3;

cr = .05*e;     % micro molar
co = 1e3;       % micro molar
k1 = 1.5e-3;        % 1 / micro molar / ms
k2 = 0.3e-3;        % 1 / ms
br = stim.BT*cr./(cr + k2/k1);
Dc = 220e-11;       % cm^2 / ms   parker
Db = 110e-11;        % cm^2 / ms    parker
F = 96485.3399;     % coulombs / mole   Faraday's constant
Kd = 1;             % micromolar
VT = 25.8;  % mV

S = spdiags([-e 2*e -e], -1:1, N, N)/dx/dx;
S(1,1) = 1/dx/dx;
S(N,N) = 1/dx/dx;
lam2 = cab.rad/2/R2;
B = lam2*S;

eN = speye(N);
R = [-Dc*S-k1*stim.BT*eN  k2*eN; k1*stim.BT*eN -Db*S-k2*eN];
Rp = 2*speye(2*N) + cab.dt*R;
Rm = 2*speye(2*N) - cab.dt*R;
[LRm,URm] = lu(Rm);

ci = cr;
Vr = fsolve(@(V) Iss(V,E,g,B,ci,co,VT),-70*e);  % initial conditions

n1 = an(Vr)./(an(Vr)+bn(Vr)); n14 = n1.^4;
m1 = am(Vr)./(am(Vr)+bm(Vr)); 
h1 = ah(Vr)./(ah(Vr)+bh(Vr)); m13h1 = m1.^3.*h1;

m1T = amT(Vr)./(amT(Vr)+bmT(Vr)); 
h1T = ahT(Vr)./(ahT(Vr)+bhT(Vr)); m12h1T = m1T.^2.*h1T;

m1N = amN(Vr)./(amN(Vr)+bmN(Vr)); 
h1N = ahN(Vr)./(ahN(Vr)+bhN(Vr)); m12h1N = m1N.^2.*h1N;

m1L = amL(Vr)./(amL(Vr)+bmL(Vr)); m12L = m1L.^2;

inaca = stim.knaca*(Na3*exp(Vr/2/VT)-(ci/co).*exp(-Vr/2/VT));

m1KCa = amKCa(Vr,cr,VT)./(amKCa(Vr,cr,VT)+bmKCa(Vr,cr,VT));

ICaT(1) = g.CaT*m12h1T(Nmid)*Phi(cr(Nmid),co,Vr(Nmid),VT);
ICaN(1) = g.CaN*m12h1N(Nmid)*Phi(cr(Nmid),co,Vr(Nmid),VT);
ICaL(1) = g.CaL*m12L(Nmid)*Phi(cr(Nmid),co,Vr(Nmid),VT);
INaCa(1) = inaca(Nmid);
ICa(1) = ICaT(1)+ICaN(1)+ICaL(1)+INaCa(1);
IKCa(1) = g.KCa*m1KCa(Nmid)*(Vr(Nmid)-E.K);
IK(1) = g.K*n14(Nmid)*(Vr(Nmid)-E.K);
ca(1) = cr(Nmid);
bca(1) = br(Nmid);
vmid(1) = Vr(Nmid);

u = [cr; br];           % initially at rest

ICa1 = (g.CaT*m12h1T + g.CaN*m12h1N + g.CaL*m12L).*Phi(u(1:N),co,Vr,VT);

dt = cab.dt;

v = Vr;
a = an(v);  b = bn(v);
n2 = ( (2/dt-a-b).*n1 + 2*a) ./ (2/dt + a + b); n24 = n2.^4;

a = am(v);  b = bm(v);
m2 = ( (2/dt-a-b).*m1 + 2*a) ./ (2/dt + a + b);

a = ah(v);  b = bh(v);
h2 = ( (2/dt-a-b).*h1 + 2*a) ./ (2/dt + a + b); m23h2 = m2.^3.*h2;

a = amT(v);  b = bmT(v);
m2T = ( (2/dt-a-b).*m1T + 2*a) ./ (2/dt + a + b);
a = ahT(v);  b = bhT(v);
h2T = ( (2/dt-a-b).*h1T + 2*a) ./ (2/dt + a + b); m22h2T = m2T.^2.*h2T;

a = amN(v);  b = bmN(v);
m2N = ( (2/dt-a-b).*m1N + 2*a) ./ (2/dt + a + b);
a = ahN(v);  b = bhN(v);
h2N = ( (2/dt-a-b).*h1N + 2*a) ./ (2/dt + a + b); m22h2N = m2N.^2.*h2N;

a = amL(v);  b = bmL(v);
m2L = ( (2/dt-a-b).*m1L + 2*a) ./ (2/dt + a + b); m22L = m2L.^2;

a = amKCa(v,cr,VT);  b = bmKCa(v,cr,VT);
m2KCa = ( (2/dt-a-b).*m1KCa + 2*a) ./ (2/dt + a + b);

inaca = stim.knaca*(Na3*exp(v/2/VT)-(u(1:N)/co).*exp(-v/2/VT));

ICaT(2) = g.CaT*m22h2T(Nmid)*Phi(u(Nmid),co,Vr(Nmid),VT);
ICaN(2) = g.CaN*m22h2N(Nmid)*Phi(u(Nmid),co,Vr(Nmid),VT);
ICaL(2) = g.CaL*m22L(Nmid)*Phi(u(Nmid),co,Vr(Nmid),VT);
IKCa(2) = g.KCa*m2KCa(Nmid)*(Vr(Nmid)-E.K);
IK(2) = g.K*n24(Nmid)*(Vr(Nmid)-E.K);
INaCa(2) = inaca(Nmid);
ICa(2) = ICaT(2)+ICaN(2)+ICaL(2)+INaCa(2);
ca(2) = u(Nmid);
bca(2) = u(Nmid+N);
vmid(2) = Vr(Nmid);

t = 0;
if pinc 
   figure(1)
   plot3(x,t*ones(N,1),u(N+1:2*N),'color','k')
   hold on
   figure(2)
   plot3(x,t*ones(N,1),u(1:N),'color','k')
   hold on
   figure(3)
   plot3(x,t*ones(N,1),-ICa1,'color','k')
   hold on
   figure(4)
   plot3(x,t*ones(N,1),v,'color','k')
   hold on
end

i1 = 0;
t = cab.dt;
i2 = Iapp*(mod(t,stim.per)<stim.pw);

dB = diag(B);

tf = 2*Cm/cab.dt;

d2 = g.Na*m23h2 + g.K*n24 + g.Cl + g.KCa*m2KCa;

ICa2 = (g.CaT*m22h2T + g.CaN*m22h2N + g.CaL*m22L).*Phi(u(1:N),co,v,VT);
   
f = g.Na*(m23h2 + m13h1)*E.Na + g.K*(n24+n14)*E.K + 2*g.Cl*E.Cl ...
      - 2*inaca + i2 + i1 - (ICa1 + ICa2) + g.KCa*(m2KCa+m1KCa)*E.K;

r = (tf - d2).*v - B*v + f;

B(1:N+1:end) = dB + d2 + tf;         % update the diagonal

v = B\r;   % v at dt

if pinc
   figure(1)
   plot3(x,t*ones(N,1),u(N+1:2*N),'color','k')
   figure(2)
   plot3(x,t*ones(N,1),u(1:N),'color','k')
   figure(3)
   plot3(x,t*ones(N,1),-ICa2,'color','k')
   figure(4)
   plot3(x,t*ones(N,1),v,'color','k')
end

i1 = i2; n1 = n2; m1 = m2; h1 = h2; d1 = d2; m13h1 = m23h2; n14 = n24; 
h1N = h2N; m1N = m2N; h1T = h2T; m1T = m2T; m1L = m2L; ICa1 = ICa2;
m1KCa = m2KCa;

for j=2:Nt,
    
    rad = cab.rad;
    
    p1 = -(2/rad)*(ICa2-inaca)/(2*F);
    p2 = (2/rad)*(stim.Jpmca*u(1:N)./(Kd+u(1:N)));
    p = p1 - p2;
    tmp = k1*u(1:N).*u(N+1:2*N);
    q = 2*[tmp + p; -tmp];
    u = URm \ ( LRm \ (Rp*u + dt*q) );
    
    t = j*dt;

    i2 = Iapp*(mod(t,stim.per)<stim.pw);

    a = an(v);  b = bn(v);
    n2 = ( (2/dt-a-b).*n1 + 2*a) ./ (2/dt + a + b); n24 = n2.^4;

    a = am(v);  b = bm(v);
    m2 = ( (2/dt-a-b).*m1 + 2*a) ./ (2/dt + a + b);    

    a = ah(v);  b = bh(v);
    h2 = ( (2/dt-a-b).*h1 + 2*a) ./ (2/dt + a + b); m23h2 = m2.^3.*h2;

    a = amT(v);  b = bmT(v);
    m2T = ( (2/dt-a-b).*m1T + 2*a) ./ (2/dt + a + b);
    a = ahT(v);  b = bhT(v);
    h2T = ( (2/dt-a-b).*h1T + 2*a) ./ (2/dt + a + b); m22h2T = m2T.^2.*h2T;

    a = amN(v);  b = bmN(v);
    m2N = ( (2/dt-a-b).*m1N + 2*a) ./ (2/dt + a + b);
    a = ahN(v);  b = bhN(v);
    h2N = ( (2/dt-a-b).*h1N + 2*a) ./ (2/dt + a + b); m22h2N = m2N.^2.*h2N;

    a = amL(v);  b = bmL(v);
    m2L = ( (2/dt-a-b).*m1L + 2*a) ./ (2/dt + a + b); m22L = m2L.^2;

    a = amKCa(v,u(1:N),VT);  b = bmKCa(v,u(1:N),VT);
    m2KCa = ( (2/dt-a-b).*m1KCa + 2*a) ./ (2/dt + a + b);

    d2 = g.Na*m23h2 + g.K*n24 + g.Cl + g.KCa*m2KCa;

    ICa2 = (g.CaT*m22h2T + g.CaN*m22h2N + g.CaL*m22L).*Phi(u(1:N),co,v,VT);
    
    inaca = stim.knaca*(Na3*exp(v/2/VT)-(u(1:N)/co).*exp(-v/2/VT));

    f = g.Na*(m23h2 + m13h1)*E.Na + g.K*(n24+n14)*E.K + 2*g.Cl*E.Cl ...
        - 2*inaca + i2 + i1 - (ICa1 + ICa2) + g.KCa*(m2KCa+m1KCa)*E.K;

    r = (2*tf - d2 + d1).*v - r + f;

    B(1:N+1:end) = dB + d2 + tf;         % update the diagonal

    v = B\r;

    if mod(j,pinc) == 0
        figure(1)
        plot3(x,t*ones(N,1),u(N+1:2*N),'color','k')
        figure(2)
        plot3(x,t*ones(N,1),u(1:N),'color','k')
        figure(3)
        plot3(x,t*ones(N,1),-ICa2,'color','k')
        figure(4)
        plot3(x,t*ones(N,1),v,'color','k')
    end

    ICaT(j+1) = g.CaT*m22h2T(Nmid)*Phi(u(Nmid),co,v(Nmid),VT);
    ICaN(j+1) = g.CaN*m22h2N(Nmid)*Phi(u(Nmid),co,v(Nmid),VT);
    ICaL(j+1) = g.CaL*m22L(Nmid)*Phi(u(Nmid),co,v(Nmid),VT);
    IKCa(j+1) = g.KCa*m2KCa(Nmid)*(v(Nmid)-E.K);
    IK(j+1) = g.K*n24(Nmid)*(v(Nmid)-E.K);
    INaCa(j+1) = inaca(Nmid);
    vmid(j+1) = v(Nmid);
    ICa(j+1) = ICaT(j+1)+ICaN(j+1)+ICaL(j+1)+INaCa(j+1);
    ca(j+1) = u(Nmid);
    bca(j+1) = u(Nmid+N);

    i1 = i2; n1 = n2; m1 = m2; h1 = h2; d1 = d2; m13h1 = m23h2; n14 = n24;
    h1N = h2N; m1N = m2N; h1T = h2T; m1T = m2T; m1L = m2L; ICa1 = ICa2;
    m1KCa = m2KCa;

end

if pinc > 0
    figure(1)
   xlabel('x  (cm)','fontsize',14,'color','k')
   ylabel('t  (ms)','fontsize',14,'color','k')
   zlabel('b  (\muM)','fontsize',14,'color','k')
   hold off
   figure(2)
   xlabel('x  (cm)','fontsize',14,'color','k')
   ylabel('t  (ms)','fontsize',14,'color','k')
   zlabel('c  (\muM)','fontsize',14,'color','k')
   hold off
   figure(3)
   xlabel('x  (cm)','fontsize',14,'color','k')
   ylabel('t  (ms)','fontsize',14,'color','k')
   zlabel('-I_{Ca}  (\muA/cm^2)','fontsize',14,'color','k')
   hold off
   figure(4)
   xlabel('x  (cm)','fontsize',14,'color','k')
   ylabel('t  (ms)','fontsize',14,'color','k')
   zlabel('v  (mV)','fontsize',14,'color','k')
   hold off
end

t = linspace(0,stim.Tfin,length(ICaT))';

return

function val = Iss(V,E,g,B,ci,co,VT)
val = B*V + g.Na*(am(V)./(am(V)+bm(V))).^3.*...
            (ah(V)./(ah(V)+bh(V))).*(V-E.Na) + ...
            g.K*(an(V)./(an(V)+bn(V))).^4.*(V-E.K) + g.Cl.*(V-E.Cl) + ...
   (g.CaT*(amT(V)./(amT(V)+bmT(V))).^2.*(ahT(V)./(ahT(V)+bhT(V))) + ...
    g.CaN*(amN(V)./(amN(V)+bmN(V))).^2.*(ahN(V)./(ahN(V)+bhN(V))) + ...
    g.CaL*(amL(V)./(amL(V)+bmL(V))).^2) .*Phi(ci,co,V,VT);

function val = an(v)
val = .01*(10-(v+71))./(exp(1-(v+71)/10)-1);

function val = bn(v)
val = .125*exp(-(v+71)/80);

function val = am(v)
val = .1*(25-(v+71))./(exp(2.5-(v+71)/10)-1);

function val = bm(v)
val = 4*exp(-(v+71)/18);

function val = ah(v)
val = 0.07*exp(-(v+71)/20);

function val = bh(v)
val = 1./(exp(3-(v+71)/10)+1);

function val = Phi(ci,co,v,VT)
e = exp(2*v/VT);
val = v.*(1-(ci/co).*e)./(1-e);

function val = amL(V)
a_L = 15.69; b_L = 81.5; c_L = 0.29; d_L = 10.86;
val = a_L*(b_L-V)./(exp((b_L-V)/10)-1);

function val = bmL(V)
a_L = 15.69; b_L = 81.5; c_L = 0.29; d_L = 10.86;
val = c_L*exp(-V/d_L);

function val = amN(V)
a_N = 0.19; b_N = 19.88; c_N = 0.046; d_N = 20.73;
e_N = 1.6e-4; f_N=48.46; g_N=39;
val = a_N*(b_N-V)./(exp((b_N-V)/10)-1);

function val = bmN(V)
a_N = 0.19; b_N = 19.88; c_N = 0.046; d_N = 20.73;
e_N = 1.6e-4; f_N=48.46; g_N=39;
val = c_N*exp(-V/d_N);

function val = ahN(V)
a_N = 0.19; b_N = 19.88; c_N = 0.046; d_N = 20.73;
e_N = 1.6e-4; f_N=48.46; g_N=39;
val = e_N*exp(-V/f_N);

function val = bhN(V)
a_N = 0.19; b_N = 19.88; c_N = 0.046; d_N = 20.73;
e_N = 1.6e-4; f_N=48.46; g_N=39;
val  = 1./(1+exp((g_N-V)/10));

function val = amT(V)
a_T = 0.2; b_T = 19.26; c_T = 0.009; d_T = 22.03;
e_T = 1e-6; f_T=16.26; g_T=29.79;
val = a_T*(b_T-V)./(exp((b_T-V)/10)-1);

function val = bmT(V)
a_T = 0.2; b_T = 19.26; c_T = 0.009; d_T = 22.03;
e_T = 1e-6; f_T=16.26; g_T=29.79;
val = c_T*exp(-V/d_T);

function val = ahT(V)
a_T = 0.2; b_T = 19.26; c_T = 0.009; d_T = 22.03;
e_T = 1e-6; f_T=16.26; g_T=29.79;
val = e_T*exp(-V/f_T);

function val = bhT(V)
a_T = 0.2; b_T = 19.26; c_T = 0.009; d_T = 22.03;
e_T = 1e-6; f_T=16.26; g_T=29.79;
val = 1./(1+exp((g_T-V)/10));

function val = amKCa(v,c,VT)
val = 0.28*c./(c+.48*exp(-1.7*v/VT));

function val = bmKCa(v,c,VT)
val = 0.48./(1+(c/.13e-3).*exp(2*v/VT));

