% Fully connected recurrent ANN (artificial neural network)
clear; clf; hold on;

% Network parameters
nn=500; % nodes
npat=50; % patterns
dt=.01; % time step
tau=1; % time constant
tend=.1; % end simulation

% random binary pattern with elements [-1;1]
pat=2*floor(2*rand(nn, npat))-1;

% weight matrix with Hebbian learning
w=pat*pat';
%w=w./npat-diag(1); % if we want to remove self couplings

for dist0=0.1:.1:.0; % loop over initial distnace
    u=pat(:,1);
    nflips=floor(dist0*nn);
    flag=zeros(nn,1);
    iflip=1;
    while iflip<nflips;
        ii=floor(nn*rand)+1;
        if flag(ii)==0; u(ii)=-u(ii); flag(ii)=1; iflip=iflip+1; end
    end
    s=tanh(u);
    x(1)=0;
    y(1)=(1-pat(:,1)'*s/(sqrt(pat(:,1)'*pat(:,1))*sqrt(s'*s)))/2;
    irec=1;
    for t=tdt:dt:tend;
        s=tanh(u);
        u=(1-dt/tau)*u+dt/tau*(w*s);
        irec=irec+1;
        x(irec)=t;
        y(irec)=(1-pat(:,1)'*s/sqrt(pat(:,1)'*pat(:,1))*sqrt(s'*s)))/2;
    end % time
    plot(x,y);
end % dist0
axis([0 tend 0 1]); xlabel('Time [tau]'); tlabel('Distance');
