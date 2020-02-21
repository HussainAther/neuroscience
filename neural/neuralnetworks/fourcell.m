% Simulate a net of four E cells and Hebbian Plasticity 

function fourcell(dt,Tfin,P)

close all

Nt = ceil(Tfin/dt);
w41 = zeros(Nt,1);  w43 = w41;		% weights to be followed

tauE = 2; 

gL = 0.3; VL = -68; Cm = 1; Wmax = 1;

WEE = [0 0 0 0;
       0.75 0 0 0;
       0 0.75 0 0;
       0.75 0 0.7 0];

w41(1) = WEE(4,1);  w43(1) = WEE(4,3);

tauP = 10; tauD =  10;
AP =  0.3; AD = 0.3;

pre = cell(4,1); post = cell(4,1);

for i=1:4,
    pre{i} = find(WEE(i,:));
    post{i} = find(WEE(:,i));
end

aE = (2*tauE-dt)/(2*tauE+dt);  bE = 2/(2*tauE+dt);

tref = 3; Vthr = -50; Vres = -70;

gE = zeros(4,1); VE = VL*ones(4,1); t = zeros(Nt,1);

spE = zeros(4,1); TE = -500*ones(4,1);

a = 2*Cm/dt; b = 2*gL*VL; c = a + gL;

for j=1:Nt-1,

    t = j*dt;

    spin = (t/P==round(t/P));

    gEo = gE;
    gE = aE*gE + bE*(WEE*spE + spin*[1 0 0 0]');

    VE = ( (a - (gL + gEo)).*VE + b) ./  (c + gE);

    refE = find( t-TE < tref );
    VE(refE) = Vres;

    spE = (VE>Vthr);

    if sum(spE)>0

        fspE = find(spE>0);
        TE(fspE) = t;
        VE(fspE) = Vres;
        subplot(2,1,1)
        plot(t,fspE,'k+')
        hold on

        nEs = length(fspE);
    
        for n = 1:nEs,
            
            i = fspE(n);     % cell i has spiked

            ipre = pre{i};

            if ~isempty(ipre)    % potentiate presyn weights

                Deltat = TE(ipre) - t;
                dWP = AP*exp(Deltat/tauP)';
                WEE(i,ipre) = dWP + WEE(i,ipre).*(1 - dWP);

            end

            ipost = post{i};

            if isempty(ipost)
               continue
            end

            Deltat = TE(ipost)-t;	% depress postsyn weights
            dWD = AD*exp(Deltat/tauD);
            WEE(ipost,i) = WEE(ipost,i).*(Wmax - dWD);
            
        end

    end % if

    w41(j+1) = WEE(4,1);  w43(j+1) = WEE(4,3);

end % for j

hold off
xlim([0 Tfin])
ylim([.5 4.5])
set(gca,'xticklabel',[])
set(gca,'tickdir','out')
box off
ylabel('cell','fontsize',14)

subplot(2,1,2)
t = linspace(0,Tfin,Nt);
plot(t,w41,'k',t,w43,'r')
box off
set(gca,'tickdir','out')
legend('W_{4,1}','W_{4,3}','location','SW')
legend(gca,'boxoff')
xlabel('t  (ms)','fontsize',14)

