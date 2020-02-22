%  A model of the NeuroVascularUnit following
%  K Dormanns, EMJ van Disseldorp, RG Brown and T David
%  Neurovascular coupling and the influence of luminal agonists
%  via the endothelium, J Theoretical Biology 364 (2015) 49-70
function [t,kstim,kp,Vm,Jcam,cam,rad] = nvu(dt,tfin,Jip3lum)

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
K = [0 0.5 0.4 0.1 0.5 0 0.1];    % Hai-Murphy Crossbridge coeeficients
stimon = 100;
stimoff = 210;

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

kp = zeros(jfin,1);     % [K^+]_perivascular space

% set resting values

ks(1) = 3; %3;    % mM, pressure is balanced, w'_a = 0
ka(1) = 134; % 115; % 103;
ns(1) = 146; % 117; %120;
na(1) = 15; % 5; %20;

% now determine pump rate, p, from Na and K flux

ENa = RToF*log(ns(1)/na(1));
EK = RToF*log(ks(1)/ka(1));
Kni = 10;    % mM
Kni32 = Kni^(3/2);
Kko = 1.5;   % mM
JNKnop = (na(1)^(3/2)/(na(1)^(3/2)+Kni32))*ks(1)/(ks(1)+Kko);
p = g.K*g.Na*(ENa-EK)/F/(2*g.Na+3*g.K)/JNKnop;
JNK = p*JNKnop;
Vabar = EK + 2*F*JNK/g.K;

Va(1) = Vabar;
%Vanow = Va(1)

% now determine Eh and ECl from

Eh = (Va(1)+ENa)/2;
ECl = Va(1);

% now finally determine remaining concentrations from 
% charge balance

Xa = 2e-5;    % micro mol
wa = 2e-7; % 1e-8;    % cm^3
Aa = 4e-2; % 1e-4;    % cm^2
ua(1) = wa/Aa;
%uanow = ua(1)

mat = [exp(-ECl/RToF) exp(-Eh/RToF); 1 1];
rhs = [ns(1)+ks(1); na(1)+ka(1)-Xa/wa];
cha = mat \ rhs;
ca(1) = cha(1);
ha(1) = cha(2);
cs(1) = ca(1)*exp(-ECl/RToF);
hs(1) = ha(1)*exp(-Eh/RToF);
%hanow = ha(1)
%hsnow = hs(1)

% resting levels for BK machinery

mBK = minfBK(Va(1));                   % BK gating variable
kp(1) = ka(1)*exp(Va(1)/RToF);         % mM,  where Vabar = E_BK
%kp(1) = 3;   % doubleduty
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

% initialize muscle variables

cam = zeros(jfin,1);     % [Ca^++]_muscle
Jcam = cam;               % ca current in muscle 
casr = zeros(jfin,1);     % [Ca^++]_sarcoplasmic reticulum
ip3m = zeros(jfin,1);     % [IP_3]_muscle
Vm = zeros(jfin,1);     % muscle, membrane potential
rad = zeros(jfin,1);     % radius of blood vessel
Fr = zeros(jfin,1);     % fraction of bridges

cam(1) = 0.2; % CAEC = cam(1);
casr(1) = 0.1;
ip3m(1) = 0.1;  % IP3EC = ip3m(1);
Vm(1) = -60; % -25;    % VEC = Vm(1);
rad(1) = 15;
thick = 1.5;
Jcam(1) = (1.29e-3)*(Vm(1)-100)/(1+exp(-((Vm(1)+24)/8.5)));   % (A.47)
%kp(1) = 3;        % perivascular K^+ concentration, mM

ym =  cam(1)^2/(cam(1)^2 + 0.13*exp(-(Vm(1)+27)/12));  % KCa gate

K(1) = 17*cam(1)^3;
K(6) = K(1);
HMA = [-(K(1)+K(2)+K(3)) K(4)-K(1) -K(1); 
         K(3) -(K(4)+K(5)) K(6);
         0  K(5) -(K(6)+K(7))];
HMb = [-K(1); 0; 0];

HMrest = HMA \ HMb;
Mp = HMrest(1);
AMp = HMrest(2);
AM = HMrest(3);
Fr(1) = AMp + AM;

% initialize endothelium

cae = zeros(jfin,1);     % [Ca^++]_endothelium
caer = zeros(jfin,1);     % [Ca^++]_endooplasmic reticulum
ip3e = zeros(jfin,1);     % [IP_3]_endothelium
Ve = zeros(jfin,1);     % endothelium, membrane potential

cae(1) = 0.1; % CAM = cae(1);
caer(1) = 0.1;
ip3e(1) = 0.1;  % IP3M = ip3e(1);
Ve(1) = -75; % -25;    % VM = Ve(1);

for j=1:jfin-1,
     
    t(j+1) = j*dt;
    
    % advance astrocylle variables
    
    ENa = RToF*log(ns(j)/na(j));
    EK = RToF*log(ks(j)/ka(j));
    ECl = -RToF*log(cs(j)/ca(j));
    Eh = -RToF*log(hs(j)/ha(j));
    ENBC = 2*Eh - ENa;
    
    na32 = na(j)^(3/2);
    JNK = F*p*(na32/(na32+Kni32))*ks(j)/(ks(j)+Kko);

    cot = (t(j+1)>stimon)*(t(j+1)<stimoff);      % cotransporter switch
 
    Va(j+1) = (g.Na*ENa + g.K*EK + g.Cl*ECl + g.BK*mBK*EBK + cot*g.NBC*ENBC - JNK)/...
              (Vabot + g.BK*mBK + cot*g.NBC);

    JNKCC1 = cot*g.NKCC1*(ENa + EK - 2*ECl);
    JKCC1 = cot*g.KCC1*(EK - ECl);
    JNBC = cot*g.NBC*(Va(j+1)-ENBC);
    
    et = (t(j+1)-stimon)*(t(j+1)>stimon);
    kstim(j+1) = 7e-6*(et*exp(-et) - ...
                   (1/10)*(t(j+1)>stimoff-10)*(t(j+1)<stimoff));   % potassium stimulus

    mBK = mBK + dt*( (minfBK(Va(j+1))-mBK)/tauBK(Va(j+1)) );
    EBK = RToF*log(kp(j)/ka(j));
    JBK = g.BK*mBK*(Va(j+1)-EBK);    
    
    krhs = (1/F)*(g.K*(Va(j+1)-EK) - 2*JNK - JNKCC1 - JKCC1 + JBK);
    uks = uks + dt*(krhs + kstim(j+1)/Aa);
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
    
    % compute muscle fluxes
    
    i2 = ip3m(j)^2;
    JIP3m = 0.23*i2/(1+i2);          % (A.42)   micro molar / s
    
    c2 = cam(j)^2;
    JSRup = 2.025*c2/(1+c2);        % (A.43)  micro molar / s
    
    cs2 = casr(j)^2;
    JCICRsr = 55*(cs2/(4+cs2))*(c2^2/(0.9^4+c2^2));   % (A.44) micro molar / s
    
    Jcampump = 0.24*cam(j)*( 1 + (Vm(j)+100)/250 );   % (A.45)
    
    JSRleak = 0.025*casr(j);    % (A.46)
    
    JCAm = (1.29e-3)*(Vm(j)-100)/(1+exp(-((Vm(j)+24)/8.5)));   % (A.47)
    Jcam(j+1) = JCAm;    % collect for plotting purpose
    
    JNaCam = (3.16e-3)*cam(j)*(Vm(j)+30)/(cam(j)+0.5);     % (A.48)
    
    Jstretchm = (6.1e-3)*(Vm(j)+18)/(1+exp(-7.4e-3*(30*rad(j)/thick-500)));  %(A.49)
    
    JNaKm = 4.32e-2;    % (A.50)
    
    JClm = (1.34e-3)*(Vm(j)+25);   % (A.51)
    
    JKCam = (4.46e-3)*ym*(Vm(j)+94);    % (A.52)
    
    Kactm =  c2/(c2 + 0.13*exp(-(Vm(j)+27)/12));    % (A.70) + table below it
    
    %tkp = kp + stim.kp*(t(j)>60)*(t(j)<120);
    EKIR = 4.5*kp(j) - 112;                    % Table 1, p 53 and (2.7)
    gKIR = exp((0.42)*kp(j) - (7.4e-2)*Vm(j) - 12.6);    % Table 1 and (2.8)
    JKIRm = (750/1970)*gKIR*(Vm(j)-EKIR);     % (A.53)
    
    JIP3deg = 0.1*ip3m(j);
    
    % advance perisynaptic potassium, careful with dimensions 
    %  going fom the mili-world to the micro-world
    
    kp(j+1) = kp(j) + dt*(JBK/F/ua(j) + JKIRm/1000)/1e-3; % (A.15) microM!
    
    % advance muscle processes
    
    ym = ym + dt*45*(Kactm-ym);      % (A.35) gating of muscle KCa channel
    
    %CAECt = CAEC + stim.caec*(t(j)>60)*(t(j)<120);
    cam(j+1) = cam(j) + dt*( JIP3m - JSRup + JCICRsr - JCAm + ...
        JNaCam - Jcampump + 0.1*Jstretchm + JSRleak - 0.05*(cam(j)-cae(j))); 
      
    casr(j+1) = casr(j) + dt*( JSRup - JCICRsr - JSRleak );
    
    %IP3ECt = IP3EC + stim.ip3ec*(t(j)>60)*(t(j)<120);
    ip3m(j+1) = ip3m(j) + dt*( -0.05*(ip3m(j)-ip3e(j)) - JIP3deg );
    
    Vm(j+1) = Vm(j) + dt*( 1970*(-JNaKm - JClm - 2*JCAm - ...
           JNaCam - JKCam - JKIRm - Jstretchm)- 0.5*(Vm(j)-Ve(j)));
    
    K(1) = 17*cam(j)^3;
    K(6) = K(1);
    Mp = Mp + dt*((K(4)-K(1))*AMp + K(1) - K(1)*AM - (K(1)+K(2)+K(3))*Mp);
    AMp = AMp + dt*(K(3)*Mp + K(6)*AM - (K(4)+K(5))*AMp);
    AM = AM + dt*(K(5)*AMp - (K(6)+K(7))*AM);
    Fr(j+1) = AMp + AM;
    
    thick = -rad(j) + sqrt(rad(j)^2 + 2*(20)*(3) + (3)^2);  % (A.78)
    EFr = 66 + Fr(j)*(233-66);
    R0 = 20 + Fr(j)*(0.6-1)*20;
    rad(j+1) = rad(j) + dt*(20/10)*(rad(j)*4/thick - EFr*(rad(j)-R0)/R0);
    
    % compute endothelial fluxes
    
    i2 = ip3e(j)^2;
    JIP3e = 0.23*i2/(1+i2);          % (A.55)   micro molar / s
    
    c2 = cae(j)^2;
    JERup = 0.5*c2/(1+c2);        % (A.56)  
    
    ce2 = caer(j)^2;
    JCICRer = 5*(ce2/(4+ce2))*(c2^2/(0.9^4+c2^2));   % (A.57) 
    
    Jcaepump = 0.24*cae(j);   % (A.58)
    
    Jstretche = (6.1e-3)*(Ve(j)+18)/(1+exp(-7.4e-3*(30*rad(j)/thick-500)));  %(A.59)
    
    JERleak = 0.025*caer(j);    % (A.60)
    
    logcae = log10(cae(j));
    Jcate = 6.6e-4*(50-Ve(j))*(1+tanh((logcae+0.18)/0.37))/2;  % (A.61)
    
    top = (logcae+0.4)*(Ve(j)+80.8) - 53.3;
    bot = 1.32e-3*(Ve(j)+53.3*(logcae+0.4)+80.8)^2 + 0.3;
    JBKe = 0.2*(1+tanh(top/bot));      % (A.63)
    
    JSKe = 0.3*(1+tanh((logcae+0.28)/0.389));    % (A.64)
    
    JKe = 6927*(Ve(j)+80)*(JBKe+JSKe);       % (A.62)
    
    JRe = 955*(Ve(j)+31.1);                 % (A.65)
    
    JIP3dege = 0.1*ip3e(j);                 % (A.66)
    
    J0e = 0.029;  % from nvu code
    
    % advance endothelial processes
    
    cae(j+1) = cae(j) + dt*( JIP3e - JERup + JCICRer - Jcaepump + ...
         JERleak + Jcate + J0e + Jstretche  - 0.05*(cae(j)-cam(j))); 
      
    caer(j+1) = caer(j) + dt*( JERup - JCICRer - JERleak );
    
    ip3e(j+1) = ip3e(j) + dt*( Jip3lum - JIP3dege - 0.05*(ip3e(j)-ip3m(j)) );
    
    Ve(j+1) = Ve(j) + dt*( -(JKe + JRe)/25.8 - 0.5*(Ve(j)-Vm(j))); 
    
end


% figure
% plot(t,ks,'k','linewidth',1.5)
% hold on
% plot(t,kp,'r','linewidth',1.5)
% plot(t,kstim/7e-7,'k','linewidth',1)
% xlabel('t  (s)','fontsize',14)
% ylabel('(mM)','fontsize',14)
% legend('[K^+]_s','[K^+]_p','k_{stim}','location','best')
% text(20,13,'(A)','fontsize',20)
% box off
% hold off
% 
% figure
% [ax,~,~] = plotyy(t,Vm,t,cam);
% xlabel('t   (s)','fontsize',14)
% ylabel(ax(1),'V_m   (mV)','fontsize',14)
% ylabel(ax(2),'[Ca^{2+}]_m   (\muM)','fontsize',14)
% 
% 
% figure
% plot(t,rad,'k');
% xlabel('t   (s)','fontsize',14)
% ylabel('rad   (\mum)','fontsize',14)

% figure(2)
% [ax,h1,h2] = plotyy(t,ua*Aa*1e7,t,Va);
% set(h1,'linewidth',1.5,'color','r')
% set(h2,'linewidth',1.5,'color','k')
% xlabel('t  (s)','fontsize',14)
% ylabel(ax(1), 'volume, w  (cm^3x1e-7)','fontsize',14,'color','r')
% ylabel(ax(2), 'V_a  (mV)','fontsize',14,'color','k')
% text(ax(2),80,-70,'(B)','fontsize',20)
% % box off

return

function val = minfBK(v)
val = (1+tanh((v+22)/14.5))/2;

function val = tauBK(v)
val = 1/(2.664*cosh((v+22)/(2*14.5)));

dt = 0.005;
tfin = 300;
Jip3lum = 0.18;
[t,kstim,kplo,Vmlo,Jcamlo,camlo,radlo] = nvu(dt,tfin,Jip3lum);
Jip3lum = 0.4;
[t,kstim,kphi,Vmhi,Jcamhi,camhi,radhi] = nvu(dt,tfin,Jip3lum);

figure(1)
plot(t,kstim*1e5/(4e-2),'k','linewidth',1.5)   % same scaling as astrocyte2
xlabel('t  (s)','fontsize',14)
ylabel('J_{stim}   (10^{-5} cm mM/s)','fontsize',14)
text(10,6,'(A)','fontsize',20)
box off

figure(2)
plot(t,kplo,'k','linewidth',1.5)
hold on
plot(t,kphi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('[K^+]_p   (mM)','fontsize',14)
text(10,11,'(B)','fontsize',20)
hold off
box off

figure(3)
plot(t,Vmlo,'k','linewidth',1.5)
hold on
plot(t,Vmhi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('V_m   (mV)','fontsize',14)
text(50,-20,'(C)','fontsize',20)
hold off
box off

figure(4)
plot(t,Jcamlo,'k','linewidth',1.5)
hold on
plot(t,Jcamhi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('J_{Ca,m}   (\mu M/s)','fontsize',14)
text(50,-0.05,'(D)','fontsize',20)
hold off
box off

figure(5)
plot(t,camlo,'k','linewidth',1.5)
hold on
plot(t,camhi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('[Ca^{2+}]_m   (\mu M)','fontsize',14)
text(50,0.8,'(E)','fontsize',20)
hold off
box off

figure(6)
plot(t,radlo,'k','linewidth',1.5)
hold on
plot(t,radhi,'r','linewidth',1.5)
legend('J_{lum}=0.18','J_{lum}=0.4','location','best')
xlabel('t  (s)','fontsize',14)
ylabel('\rho   (\mu m)','fontsize',14)
text(250,26,'(F)','fontsize',20)
hold off
box off
