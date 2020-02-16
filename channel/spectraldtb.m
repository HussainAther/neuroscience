function D =  spectral_dtb(drift,t,Bup,Blo,y,y0,notabs_flag)
% D =  spectral_dtb(drift,t,Bup,Blo,y,y0)
% Spectral solutions to bounded drift diffusion (drift to bound).
% This can handle arbitrary changing bounds.
%
% Inputs (all in SI units)
% ~~~~~~~~~~~~~~~
% drift:     vector of drift rates
%            a second column can be used to specify coherence so that we
%            can calculaye the probability of a correct choice - if not 2nd
%            column then uses the sign of the drift to determin D.corr
% t:         time series in seconds
% Bup & Blo: vector bounds with  Blo(t) < Bup(t)
%            or can be scalars if bounds are flat,
%            +/-Inf bounds are allowed
% y:         vector of values of y to propagate: length must be a power of 2
%            sensible range to choose is (with dt=sampling interval)
%            [min(Blo)-max(drift)*dt-4*sqrt(dt)  max(Bup)+max(drift)*dt+4*sqrt(dt)]  
% y0:        vector of initial pdf (need not sum to 1)
% notabs_flag: flag as to whether (1) or not (0) to calculate the notabs pdf (can take a
%            lot of memory) - default is 0 if  not specified 
%
% Outputs
% ~~~~~~~~~~~~~~~
% Returns D, a structure - the first four have "lo" vesions too
% D.up.p(drift)       total probability of hitting the upper bound for each drift level
% D.mean_t           mean decision time 
% D.up.mean_t(drift)  mean decision time for upper bound
% D.up.pdf_t(t,drift) probability of upper bound hit at each time (sums to D.up.p)
% D.up.cdf_t(t,drift) cumulative probability of upper bound hit at each time (ends at D.up.p)
% D.corr(drift)       probability of correct response - uses drift or 2nd
%                     column
% D.drifts            returns the drifts used
% D.bounds            returns the bounds used
% D.t                 returns the times used
%
% D.notabs.pdf(drifts,y,t) probability of not being absorbed and being at y at time t
% D.notabs.pos_t(drifts,t) probability of not being absorbed and being at y>0
% D.notabs.neg_t(drifts,t) probability of not being absorbed and being at y<0
% D.notabs.y the sets of y's considered for the pdf

if nargin<7
    notabs_flag=0; % flag to detetrmine whether to store the notabs pdf
end

if isvector(drift)
    drift_dir=sign(drift);
else
    drift_dir=sign(drift(:,2));
    drift=drift(:,1);
end

nt=length(t);
dt=t(2)-t(1); %sampling interval
nd=length(drift);
ny=length(y);

if round(log2(ny))~=log2(ny)
    error('Length of y must be a power of 2');
end

if numel(y0)~=numel(y)
    error('Length of y must be sames as y0');
end

D.bounds=[Bup Blo];
D.drifts=drift;
D.y=y;
D.t=(1:nt)*dt;
D.dt=dt';

%expand any flat bounds
if numel(Bup)==1 Bup=Bup*ones(nt,1);end
if numel(Blo)==1 Blo=Blo*ones(nt,1);end

%create fft of unit variance zero mean Gaussian
kk=[0:ny/2 -ny/2+1:-1]';
sigma=1.0;
omega=2*pi*kk/range(y);
E1=exp(-0.5*dt*sigma^2*omega.^2) ;  %fft of the normal distribution - scaled suitably by dt

% this is how to create the fft non-analytcially and is useful for other
% distributions you might want to try.
% p=normpdf(y,0,sqrt(dt));
% p=p/sum(p);
% E1=real(fft(fftshift(p)));

D.up.pdf_t=zeros(nt,nd);
D.lo.pdf_t=zeros(nt,nd);

if notabs_flag
    D.notabs.pdf=zeros(nd,ny,nt);
end

p_threshold=0.001; %threshold for proportion un-terminated to stop simulation

for q=1:nd  % iterate over drifts
    %set u to initial state
    u=y0;
    
    %shift mean of gaussian by drift
    E2=E1.*exp(-1i*omega*drift(q)*dt);
    
    pup=zeros(nt,1);
    plo=zeros(nt,1);
    
    for k=1:nt %iterate over time
        %fft our current pdf
        ufft=fft(u);
        
        %convolve with gaussian with drift in the frequency domain
        %by performing pointwide multiplication
        ufft=E2.*ufft;
        
        %turn back into time domain
        u=real(ifft(ufft));
        u=max(u,0);
        
        %select density that has crossed bounds
        pup(k)=sum(u(y>=Bup(k)));
        plo(k)=sum(u(y<=Blo(k)));
        
        D.up.pdf_t(k,q)=pup(k);
        D.lo.pdf_t(k,q)=plo(k);
        
        %keep only density within bounds
        u=u.*(y>Blo(k) & y<Bup(k));
        
        if notabs_flag
            D.notabs.pdf(q,:,k)=u;
        end
        
        %exit if our threshold is reached
        if sum(u)<p_threshold, break, end
    end
    
    if notabs_flag
        D.notabs.pos_t(:,q)=sum(D.notabs.pdf(:,y>=0)')';
        D.notabs.neg_t(:,q)=sum(D.notabs.pdf(:,y<0)')';
    end
    
    D.up.p(q)=sum(pup);
    D.lo.p(q)=sum(plo);
    
    D.up.mean_t(q)=t'*pup/D.up.p(q);
    D.lo.mean_t(q)=t'*plo/D.lo.p(q);
end

D.up.cdf_t=cumsum(D.up.pdf_t);
D.lo.cdf_t=cumsum(D.lo.pdf_t);
D.corr=D.up.p.*(drift_dir>0) + D.lo.p.*(drift_dir<0)  +0.5*(drift_dir==0);
D.mean_t=D.up.p.*D.up.mean_t+D.lo.p.*D.lo.mean_t;



