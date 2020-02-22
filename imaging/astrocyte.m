%  astocyte cell from 
%  K Dormanns, EMJ van Disseldorp, RG Brown and T David
%  Neurovascular coupling and the influence of luminal agonists
%  via the endothelium, J Theoretical Biology 364 (2015) 49-70
%
%  usage:       astrocyte2(dt,tfin)
%  example:     astrocyte2(.01,100)
%
function astrocyte2(dt,tfin)

g = struct('K',1,...
           'Na',0.1,...
           'Cl',0.1,...
           'NKCC1',2e-3,...
           'NBC',8e-2,...
           'KCC1',7e-2,...
           'BK',0.1);    % mS/cm^2   

F = 9.649e4;  % Faraday's constant
RToF = 25.8;     % mV
Lp = 2e-8;    % cm/s/mM

jfin = ceil(tfin/dt);
na = zeros(jfin,1);     % [Na^+]_astrocyte
ns = zeros(jfin,1);     % [Na^+]_synaptic cleft
ka = zeros(jfin,1);     % [K^+]_astrocyte
ks = zeros(jfin,1);     % [K^+]_synaptic cleft
ha = zeros(jfin,1);     % [HCO_3^-]_astrocyte
hs = zeros(jfin,1);     % [HCO_3^-]_synaptic cleft
ca = zeros(jfin,1);     % [Cl^-]_astrocyte
cs = zeros(jfin,1);     % [Cl^-]_synaptic cleft
ua = zeros(jfin,1);     % astrocyte, volume/area
Va = zeros(jfin,1);     % astrocyte, membrane potential
t = zeros(jfin,1);      % time
kstim = t;

kp = zeros(jfin,1);     % [K^+]_perisynaptic space

% set resting values

ks(1) = 3; %3;    % mM, pressure is balanced, w'_a = 0
ka(1) = 134; % 115; % 103;
ns(1) = 146; % 117; %120;
na(1) = 15; % 5; %20;

% now determine pump rate, p, from Na and K flux

ENa = RToF*log(ns(1)/na(1))
EK = RToF*log(ks(1)/ka(1))
Kni = 10;    % mM
Kni32 = Kni^(3/2);
Kko = 1.5;   % mM
JNKnop = (na(1)^(3/2)/(na(1)^(3/2)+Kni32))*ks(1)/(ks(1)+Kko);
p = g.K*g.Na*(ENa-EK)/F/(2*g.Na+3*g.K)/JNKnop
JNK = p*JNKnop;
Vabar = EK + 2*F*JNK/g.K

Va(1) = Vabar;
%Vanow = Va(1)

% now determine Eh and ECl from

Eh = (Va(1)+ENa)/2
ECl = Va(1)

% now finally determine remaining concentrations from 
% charge balance

Xa = 2e-5;    % micro mol
wa = 2e-7; % 1e-8;    % cm^3
Aa = 4e-2; % 1e-4;    % cm^2
ua(1) = wa/Aa;
uanow = ua(1)

mat = [exp(-ECl/RToF) exp(-Eh/RToF); 1 1];
rhs = [ns(1)+ks(1); na(1)+ka(1)-Xa/wa];
cha = mat \ rhs;
ca(1) = cha(1);
ha(1) = cha(2);
cs(1) = ca(1)*exp(-ECl/RToF);
hs(1) = ha(1)*exp(-Eh/RToF);
hanow = ha(1)
hsnow = hs(1)

% resting levels for BK machinery

mBK = minfBK(Va(1));                   % BK gating variable
kp(1) = ka(1)*exp(Va(1)/RToF);         % mM,  where Vabar = E_BK
EBK = RToF*log(kp(1)/ka(1));

% prepare to march

us = ua(1)/2;      % scaled volume of synaptic cleft

uks = us*ks(1);
uka = ua(1)*ka(1);
uns = us*ns(1);
una = ua(1)*na(1);
uhs = us*hs(1);
uha = ua(1)*ha(1);

Vabot = g.Na + g.K + g.Cl;

for j=1:jfin-1,
     
    t(j+1) = j*dt;
    ENa = RToF*log(ns(j)/na(j));
    EK = RToF*log(ks(j)/ka(j));
    ECl = -RToF*log(cs(j)/ca(j));
    Eh = -RToF*log(hs(j)/ha(j));
    ENBC = 2*Eh - ENa;
    
    na32 = na(j)^(3/2);
    JNK = F*p*(na32/(na32+Kni32))*ks(j)/(ks(j)+Kko);

    cot = (t(j+1)>10)*(t(j+1)<50);      % cotransporter switch
 
    Va(j+1) = (g.Na*ENa + g.K*EK + g.Cl*ECl + g.BK*mBK*EBK + cot*g.NBC*ENBC - JNK)/...
              (Vabot + g.BK*mBK + cot*g.NBC);

    JNKCC1 = cot*g.NKCC1*(ENa + EK - 2*ECl);
    JKCC1 = cot*g.KCC1*(EK - ECl);
    JNBC = cot*g.NBC*(Va(j+1)-ENBC);
    
    et = (t(j+1)-10)*(t(j+1)>10);
    kstim(j+1) = 7e-6*(et*exp(-et) - ...
                   (1/10)*(t(j+1)>40)*(t(j+1)<50));   % potassium stimulus

    mBK = mBK + dt*( (minfBK(Va(j+1))-mBK)/tauBK(Va(j+1)) );
    EBK = RToF*log(kp(j)/ka(j));
    JBK = g.BK*mBK*(Va(j+1)-EBK);
    kp(j+1) = kp(j) + dt*(JBK/F/ua(j)/1e-3);
    
    krhs = (1/F)*(g.K*(Va(j+1)-EK) - 2*JNK - JNKCC1 - JKCC1 + JBK);
    % uks = uks + dt*(krhs + kstim(j+1)/Aa);  % this has JBK
    uks = uks + dt*(krhs-JBK/F + kstim(j+1)/Aa); % better
    uka = uka - dt*krhs;
    
    nrhs = (1/F)*(g.Na*(Va(j+1)-ENa) + 3*JNK - JNKCC1 - JNBC);
    uns = uns + dt*(nrhs - kstim(j+1)/Aa);
    una = una - dt*nrhs;
    
    uhs = uhs - dt*2*(1/F)*JNBC;
    uha = uha + dt*2*(1/F)*JNBC;
    
    ua(j+1) = ua(j) + dt*2*Lp*( na(j)-ns(j) + ka(j)-ks(j) );
    
    us = (3/2)*ua(1) - ua(j+1);
    
    ns(j+1) = uns/us;            % recover concentrations
    na(j+1) = una/ua(j+1);
    ks(j+1) = uks/us;
    ka(j+1) = uka/ua(j+1);
    hs(j+1) = uhs/us;
    ha(j+1) = uha/ua(j+1);
    
    ca(j+1) = na(j+1) + ka(j+1) - ha(j+1) - Xa/Aa/ua(j+1);  % astro charge balance
    cs(j+1) = ns(j+1) + ks(j+1) - hs(j+1);          % syn cleft charge balance
    
end

% figure(1)
% plot(t,ks,'k','linewidth',1.5)
% hold on
% plot(t,kp,'r','linewidth',1.5)
% plot(t,kstim/7e-7,'k','linewidth',1)
% xlabel('t  (s)','fontsize',14)
% ylabel('(mM)','fontsize',14)
% legend('[K^+]_s','[K^+]_p','k_{stim}','location','best')
% text(60,13,'(A)','fontsize',20)
% box off
% hold off

figure(1)
[ax,h1,h2]=plotyy(t,kp,t,1e5*kstim/Aa)
hold on
plot(t,ks,'r--','linewidth',1.5,'parent',ax(1))
%plot(t,kstim/7e-7,'k','linewidth',1)
set(ax(1),'ycolor',[1 0 0])
set(h1,'linewidth',1.5,'color','r')
set(ax(2),'ycolor',[0 0 0])
set(h2,'linewidth',1.5,'color','k')
xlabel('t  (s)','fontsize',14)
ylabel(ax(1),'(mM)','fontsize',14,'color','r')
ylabel(ax(2),'J_{stim}  (10^{-5}cm mM/s)','fontsize',14,'color','k')
legend('[K^+]_p','[K^+]_s','location','best')
text(60,13,'(A)','fontsize',20)
%box off
hold off

figure(2)
[ax,h1,h2] = plotyy(t,ua*Aa*1e7,t,Va)
set(ax(1),'ycolor',[1 0 0])
set(h1,'linewidth',1.5,'color','r')
set(ax(2),'ycolor',[0 0 0])
set(h2,'linewidth',1.5,'color','k')
xlabel('t  (s)','fontsize',14)
ylabel(ax(1), 'volume, w_a  (10^{-7} cm^3)','fontsize',14,'color','r')
ylabel(ax(2), 'V_a  (mV)','fontsize',14,'color','k')
text(60,-63,'(B)','fontsize',20,'parent',ax(2))
% box off

return

function val = minfBK(v)
val = (1+tanh((v+22)/14.5))/2;

function val = tauBK(v)
val = 1/(2.664*cosh((v+22)/(2*14.5)));
