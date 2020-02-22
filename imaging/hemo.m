% Compute the hemodynamic response following
%  Buxton et al., NeuroImage 2004
function hemo

close all
gethemo(0)   % without neural inhibition
gethemo(2)   % with neural inhibition

return

function gethemo(kappa)

if kappa==0
    figure(1)
else
    figure(2)
end

Tfin = 110;    % s
dt = 0.1;      % s

Nt = Tfin/dt;
stimE = zeros(Nt,1);
stimI = stimE;
N0 = 0;
Nresp = N0*ones(Nt,1);
%kappa = 0; %2;
tauI = 3;

for j=2:Nt,

    t = (j-1)*dt;
    if abs(t-1.5)<0.5 || abs(t-21.5)<0.5 || abs(t-41.5)<0.5 || abs(t-43.5)<0.5 || abs(t-81)<10
       stimE(j) = 1;
    end
    %Nresp(j) = max(stimE(j)-stimI(j),-N0);
    stimI(j) = stimI(j-1) + dt*(kappa*Nresp(j-1)-stimI(j-1))/tauI;
    Nresp(j) = N0 + max(stimE(j)-stimI(j),-N0);

end

tim = (dt:dt:Tfin)';
plot(tim,stimE+4.5,'k','linewidth',1.5)
hold on
if (kappa==0)
    text(5,5.8,'Excitatory Stimulus','fontsize',14)
    text(100,5.25,'(A)','fontsize',20)
else
    plot(tim,stimI+4.5,'r','linewidth',1.5)
    text(5,5.8,'Excitatory (black) and Inhibitory (red) Stimulus','fontsize',14)
    text(100,5.25,'(B)','fontsize',20)
end

plot(tim,Nresp/max(Nresp) + 3,'k','linewidth',1.5)
text(10,4.2,'Neural Response','fontsize',14)

hat = h(tim);
hatd = [zeros(10,1); hat(1:end-10)];
fhat = fft(hat);
fhatd = fft(hatd);
fNresp = fft(Nresp);
f = ifft(fhat.*fNresp);
fd = ifft(fhatd.*fNresp);
rf = real(f);
rfd = real(fd);

%plot(tim,rfd,'b:','linewidth',1.5);

% predict response to double stim

rf1 = rf(1:200);
rf2 = rf1+[zeros(20,1); rf1(1:180)];

% predict response to 20 stims

rf20 = [rf1; zeros(200,1)];
for d=1:19,
    tlo = d*10+1; % d*10+1;
    thi = 200+d*10; % 200+d*10;
    rf20(tlo:thi) = rf20(tlo:thi) + rf1;
end
plot(tim(701:1100),rf20/max(rf20) + 1.5,'k--','linewidth',1.5)
plot(tim,rf/max(rf20) + 1.5,'r','linewidth',1.5);
plot(tim(401:600),rf2/max(rf20) + 1.5,'k--','linewidth',1.5)
text(5,2,'Flow Response','fontsize',14)

f = 1 + 0.5*rf/10;
m = 1 + 0.25*rf/10;
bold = 1-m.*(f.^(-0.6));

% predict response to double stim

bold1 = bold(1:200);
bold2 = bold1+[zeros(20,1); bold1(1:180)];

% predict response to 20 stims

bold20 = [bold1; zeros(200,1)];
for d=1:19,
    tlo = d*10+1;
    thi = 200+d*10;
    bold20(tlo:thi) = bold20(tlo:thi) + bold1;
end
plot(tim(701:1100),bold20/max(bold20),'k--','linewidth',1.5)
plot(tim,bold/max(bold20),'r','linewidth',1.5)
plot(tim(401:600),bold2/max(bold20),'k--','linewidth',1.5)
text(5,0.5,'BOLD Response','fontsize',14)

hold off
box off
xlim([0 110])
xlabel('time  (s)','fontsize',14)
set(gca,'ytick',[])


function val = h(t)
val = (t.^3).*exp(-t)/6;

