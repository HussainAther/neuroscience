% Simulate a net of excitatory and inhibitory cells 
% In order to test the strength of this association we measure the learned network’s ability to complete 
% incomplete input. In particular, we systematically drop input spikes and count the average number of 
% dropped output spikes. 
%
%  usage:    EInetComp(N,dt,Tfin,P,den)
%
%  where:    N = struct('E',# of E cells,'I',# of I cells)
%           dt = timestep (ms)
%         Tfin = final time (ms)
%            P = period of stimulus (ms)
%          den = struct('EE',den of EE cnxs,'EI',lw,'II',lw,'IE') 
%
%  e.g.:   N = struct('E',80,'I',20);
%          den = struct('EE',.25,'EI',.25,'IE',.25,'II',.05);
%          EInet(N,.02,90,40,den)
%
%	EInetComp(N,.02,90,40,den,'WaNov16')
%

function spEcum = EInetComp(N,dt,Tfin,P,den,fname)

close all

Tfin = 10;
Nt = ceil(Tfin/dt);

tauE = 2; VsynE = 0;
tauI = 1; VsynI = -70;

gL = 0.3; VL = -68; Cm = 1;

if nargin == 6
   
   load(fname)
   WEE = W(1:N.E,1:N.E);
   WII = W(N.E+1:end,N.E+1:end);
   WIE = W(1:N.E,N.E+1:end);
   WEI = W(N.E+1:end,1:N.E);

else

   WEE = sprand(N.E,N.E,den.EE);
   WEE = (WEE - spdiags(diag(WEE),0,N.E,N.E))/5;

   WEI = sprand(N.I,N.E,den.EI)/4;

   WIE = sprand(N.E,N.I,den.IE)/4;

   WII = sprand(N.I,N.I,den.II);
   WII = (WII - spdiags(diag(WII),0,N.I,N.I))/4;

   W = [WEE WIE; WEI WII];

end

Ntot = N.E+N.I;

aE = (2*tauE-dt)/(2*tauE+dt); bE = 2/(2*tauE+dt);
aI = (2*tauI-dt)/(2*tauI+dt); bI = 2/(2*tauI+dt);
tref = 3; Vthr = -50; Vres = -70;
 
a = 2*Cm/dt; b = 2*gL*VL; c = a + gL;
WinEE = zeros(N.E,1); 

Nummiss = 3;
spEcum = cell(Nummiss+1,1);
nout = zeros(Nummiss+1,1);

for nummiss = 0:Nummiss,
    
    nexp = factorial(16)/(factorial(nummiss)*factorial(16-nummiss));
    %pause
    
    for nn = 1:nexp,

        spEcum{nummiss+1} = [];

        cellmiss = ceil(rand(nummiss)*16);
        WinEE(1:ceil(N.E*.2)) = 1;
        WinEE(cellmiss) = 0;
        
        gEE = zeros(N.E,1); VE = VL*ones(N.E,1); t = zeros(Nt,1);
        gII = zeros(N.I,1); VI = VL*ones(N.I,1);
        gEI = zeros(N.I,1);
        gIE = zeros(N.E,1);
        spE = zeros(N.E,1); TE = -10*ones(N.E,1);
        spI = zeros(N.I,1); TI = -10*ones(N.I,1);

        for j=1:Nt-1,

            t = j*dt;
            spin = (t==1);
            gEEo = gEE;
            gEE = aE*gEE + bE*(WEE*spE + spin*WinEE);
            gEIo = gEI;
            gEI = aE*gEI + bE*WEI*spE;
            gIEo = gIE;
            gIE = aI*gIE + bI*WIE*spI;
            gIIo = gII;
            gII = aI*gII + bI*WII*spI;
            VE = ( (a - (gL + gEEo + gIEo)).*VE + b + (gIE+gIEo)*VsynI ) ./  (c + gEE + gIE);
            VI = ( (a - (gL + gIIo + gEIo)).*VI + b + (gII+gIIo)*VsynI ) ./  (c + gII + gEI);
            refE = find( t-TE < tref );
            VE(refE) = Vres;
            spE = (VE>Vthr);
            if sum(spE)>0
                fsp = find(spE>0);
                spEcum{nummiss+1} = union(spEcum{nummiss+1},fsp);
                TE(fsp) = t;
                VE(fsp) = Vres;
                plot(t,fsp,'k.')
                hold on
            end
            refI = find( t-TI < tref );
            VI(refI) = Vres;
            spI = (VI>Vthr);
            if sum(spI)>0
                fsp = find(spI>0);
                TI(fsp) = t;
                VI(fsp) = Vres;
%                 plot(t,N.E+fsp,'r.')
%                 hold on
            end
        end  % for j
        
        hold off
       
        %axis([0 6 1 80])
        %drawnow

        spEcum{nummiss+1} = setdiff(spEcum{nummiss+1},[1:16]);
        %celldisp(spEcum)
        %pause

        outmiss = setdiff(spEcum{1},spEcum{nummiss+1});
        nout(nummiss+1) = nout(nummiss+1) + length(outmiss)/length(spEcum{1});
        
        %ylim([0 Tfin]);
        %set(gca,'ytick',[1:N])
        set(gca,'tickdir','out')
        box off
        xlabel('t  (ms)','fontsize',14)
        ylabel('cell','fontsize',14)

    end % nn
    
    nout(nummiss+1) = nout(nummiss+1)/nexp
    %pause
    
end % nummiss

