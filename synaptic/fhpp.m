%   Examine the phase plane behavior of the FitzHugh reduction
%   of HH (Hodgkin Huxley Hodgkin-Huxley)
%
%    fhpp(I,dt,Tfin,I0)      I = steady current,  I0 = perturbation
%
%    e.g.:  fhpp(20,.01,50,1)   or   fhpp(30,.01,50,1)
%

function fhpp(I,dt,Tfin,I0)

close all

A = 4*pi*(1e-6);% (cm)^2

E = struct('K', -77, 'Na', 56, 'Cl', -68); % reversal potentials, mV

G = struct('K', 36, 'Na', 120, 'Cl', 0.3);  % channel conductances, mS/cm^2

V = E.K+.1:.09:E.Na-.5;

a = -(V-E.K)*G.K;
da = -G.K;

ma = am(V); mb = bm(V);
m = ma./(ma+mb);
dm = (dam(V).*mb - ma.*dbm(V))./(ma+mb).^2;
m3 = m.^3;
dm3 = 3*m.^2.*dm;

b = G.Na*m3.*(V-E.Na);
db = G.Na*(m3 + dm3.*(V-E.Na));

%for I=0:5:25

IoA = I*1e-6/A;

c = -G.Na*m3.*(V-E.Na)*0.7 - G.Cl*(V-E.Cl) + IoA;
dc = -0.7*db - G.Cl;

ninf = an(V)./(an(V)+bn(V));

thc = da.*ninf.^4 + db.*ninf.^2 + dc - an(V) - bn(V);

Ith = -4*pi*(a.*ninf.^4 + b.*ninf.^2 + c-IoA);

zp = (-b + sqrt(b.^2-4*a.*c))/2./a;
zm = (-b - sqrt(b.^2-4*a.*c))/2./a;

figure(1)
plot(V,real(zp),'k',V,imag(zp),'k--')
hold on
plot(V,real(zm),'r',V,imag(zm),'r--')
hold off

figure(2)
n = sqrt(real(zm))
plot(V,n,'k')
hold on
plot(V,ninf,'r')
hold off
xlabel('V  (mV)','fontsize',14)
ylabel('n','fontsize',14)

figure(3)
plot(V,Ith,'k')
xlabel('V  (mV)','fontsize',14)
ylabel('I  (pA)','fontsize',14)
xlim([-75 -66])
grid

figure(5)
plot(V,thc,'k')
xlabel('V  (mV)','fontsize',14)
ylabel('S  (1/ms)','fontsize',14)


%end

Vr = fsolve(@(V) Iss(V,E,G,IoA),-71)   % find rest potential 

Nt = ceil(Tfin/dt);
t = zeros(Nt,1); V = t; n = V;
t(1) = 0;
V(1) = Vr;
n(1) = an(Vr)/(an(Vr)+bn(Vr));  

td = 1/dt;

for j = 2:Nt;

      t(j) = (j-1)*dt;

      Istim = I0*(t(j)-dt/2>2)*(1e-6);
      
      a = an(V(j-1));  b = bn(V(j-1)); c = (a+b)/2;
      n(j) = ( (td-c)*n(j-1) + a ) / (td + c);
      
      cK = G.K*n(j)^4;
      
      a = am(V(j-1));  b = bm(V(j-1)); 
      m = a / (a + b);
      
      %h = 0.87 - n(j);
      h = 0.7 - n(j)^2;
      
      cNa = G.Na*m^3*h;
      
      top = 2*V(j-1)*td + cK*E.K + cNa*E.Na + G.Cl*E.Cl + IoA + Istim/A;
      
      bot = 2*td + cK + cNa + G.Cl;
      
      Vmid = top/bot;

      V(j) = 2*Vmid - V(j-1); 

end

figure(4)
subplot(2,1,1)
plot(t,V,'k')
%ylim([-82 60])
box off
xlabel('t  (ms)','fontsize',14)
ylabel('V  (mV)','fontsize',14)
%axis tight

subplot(2,1,2)
plot(V,n,'k')
xlabel('V  (mV)','fontsize',14)
ylabel('n','fontsize',14)
box off
%axis image

function val = Iss(V,E,G,IoA)
m = am(V)/(am(V)+bm(V));
n = an(V)/(an(V)+bn(V));
%h = 0.87 - n;
h = 0.7 - n.^2;
val = G.Na*m^3*h*(V-E.Na) + G.K*n^4*(V-E.K) + G.Cl*(V-E.Cl)-IoA;

function val = an(V)
val = .01*(10-(V+71))./(exp(1-(V+71)/10)-1);

function val = bn(V)
val = .125*exp(-(V+71)/80);

function val = am(V)
val = .1*(25-(V+71))./(exp(2.5-(V+71)/10)-1);

function val = bm(V)
val = 4*exp(-(V+71)/18);

function val = ah(V)
val = 0.07*exp(-(V+71)/20);

function val = bh(V)
val = 1./(exp(3-(V+71)/10)+1);

function val = dam(v)
tmp = exp(-(46+v)/10);
val = -( tmp.*(56+v) - 10 )./(tmp-1).^2/100;

function val = dbm(v)
val = -(2/9)*exp(-(v+71)/18);

