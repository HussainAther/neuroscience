% Triple State Excitable Element Neuron Model
%   constant superthreshold input

% INPUTS   ==============================================================
   NT = 101;        % number of time steps
   pg = 0.5;        % transition probability
   uS = 4;          % membrane potential scaling 
   tStep = 10e-3;   % 1 time step = 10 ms

   
% CALCULATIONS   ========================================================
   u = zeros(NT,1);    % membrane potential
   t = 0:NT-1;         % time steps
   t = t .* tStep;     % time in seconds
   u(1) = 0;          % initial condition:  activated state

for c = 1 : NT-1
   if u(c) ==  0;   flagS = 0;  end
   if u(c) == uS;   flagS = 1;  end     
   if u(c) == -1;   flagS = -1; end
    
   switch flagS
     case 0
        u(c+1) = uS;
     case 1
        u(c+1) = -1;
     case -1    
        R = rand(1,1); 
             if R <= pg
               u(c+1) = 0;
             else
               u(c+1) = -1;
             end
   end  % end switch

end     % end for loop

% Average lieftime Tavg and transition probability pt
    pt = linspace(0.1,0.95,1000);   
    Tavg = -1 ./log(1-pt);
 
    fR = 1 ./ ((3+Tavg).*tStep);
 
% COMMAND WINDOW OUTPUT ==================================================
   lambda = -log(1-pg);      % decay constant
      
   disp('  ');
   fprintf('Transition probability from refractory to resting state  pg =  %2.2f \n',pg);
   disp('  ');
   Nactivated = length(u(u==uS));   % number of time steps for activated state
   fprintf('Number of time steps for activation state =  %2.0f \n',Nactivated);
   disp('  ');
   Nrefractory = length(u(u==-1));   % number of time steps for refractory state
   fprintf('Number of time steps for refractory state =  %2.0f \n',Nrefractory);
   disp('  ');
   Nrest = length(u(u==0));   % number of time steps for rest state
   fprintf('Number of time steps for rest state =  %2.0f \n',Nrest);
   disp('  ');
   firingRate = (Nactivated-1)/ max(t);   % firing rate of neuron
   fprintf('        Firing rate of neuron =  %2.3f Hz \n',firingRate);
   
   disp('  ');
   tAvg = 1 / lambda;   % average lifetime in refractory state
   firingPeriod = (tAvg + 3) * tStep;       % Theory: firing period [s]
   firingFrequency = 1 / firingPeriod;   % Theory: firing frequency {Hz]
   fprintf('Theory: Firing rate of neuron =  %2.3f  Hz\n',firingFrequency);
   
   disp('  ');
   fprintf('Theory: Average lifetime in refractory state =  %2.3f time steps \n',tAvg);
  
   
   disp(' new')
   [peaks, locs] = findpeaks(u,'Threshold',4);
numPeaks = length(peaks);
fPeaks1 = numPeaks/NT;

dt = zeros(numPeaks-1,1);
for c = 1:numPeaks-1
   dt(c) = locs(c+1) - locs(c);
end

fPeaks2 = 1 / mean(dt);

tAvg = -1/log(1-pg) + 3;
fPeaks3 = 1/tAvg;

fprintf('spike rate #1  =  %2.3f \n',fPeaks1);
disp('  ')
fprintf('spike rate #2  =  %2.3f \n',fPeaks2);
disp('  ')
fprintf('spike rate #3  =  %2.3f \n',fPeaks3);

       
% GRAPHICS ===============================================================

figure(1)
set(gcf,'units','normalized','Position',[0.1 0.32 0.22,0.25]);
set(gca,'fontsize',14);
plot(t,u,'lineWidth',2);
xlabel('time  t  [s]');
ylabel('membrane potential  u  [a.u.]');
axis([0 max(t) -1.2 uS+0.2]);
grid on
tm1 = 'p_g  =  ';
tm2 = num2str(pg, '%2.2f \n');
tm = [tm1 tm2];
title(tm);

figure(2)
set(gcf,'units','normalized','Position',[0.4 0.32 0.22,0.25]);
set(gca,'fontsize',14);
plot(pt,Tavg,'lineWidth',2);
xlabel('transition probability  p_g');
ylabel('average lifetime  t_{avg}  [time steps]');
%axis([0 max(t) -1.2 uS+0.2]);
grid on

figure(3)
set(gcf,'units','normalized','Position',[0.7 0.32 0.22,0.25]);
set(gca,'fontsize',14);
plot(pt,fR,'lineWidth',2);
xlabel('transition probability  p_g');
ylabel('firing rate   f_R   [Hz]');
%axis([0 max(t) -1.2 uS+0.2]);
grid on


