%{
For a Fitzhugh-Nagumo traveling wave, let the number of neurons be given by N. A naı¨ve approach might be to
allow the initial conditions to be a 2N matrix with the first row representing the N initial
voltages and the second row the N values of the initial recovery variables. However, recall
that the ode45 solver will accept only a single column array for the initial conditions. What
you will do to satisfy this requirement is let the initial value column vector be of length 2N
and let the first N elements represent the initial voltages of the N neurons and the last N
elements represent the initial recovery values. The ODE solver will produce as its output
a t1 time vector, and a t2N matrix, whose first N columns represent the evolution of
the voltage of the population of neurons as time progresses and whose second N columns
represent the evolution of the recovery variables.
The stimulus to initiate the wave that you will use is I ¼ 6 for the first 0.5s; then the stimulus 
will be off, I ¼ 0, for the rest of the time. You can choose where along the line of
neurons to initiate the wave. In the following example, you stimulate the center cells.
The following script, FNmain.m, will produce a traveling wave of activity along a onedimensional population of N neurons whose dynamics are governed by the FitzhughNagumo equations, as shown in Figure 23.3.
}%

global N I BC %by making these variables global they can exist within
%the workspace of functions without explicitly being input to the functions
N=128; %number of neurons
v0(1:N)=-1.5; %initial conditions for V variable
v0(N+1:2*N)=-3/8; %initial conditions for R variable
I=6; %the input stimulus value
BC = 2; %set to 1 (free) or 2 (periodic boundary conditions)
%
tspan=[0:.1:.5]; %time with stimulus
[t1,v1] = ode45('F_N',tspan,v0);
tspan=[.5:.1:25]; % time without stimulus
I=0;%turn off stimulus
[t2,v2] = ode45('F_N',tspan,v1(end,:)'); %note: initial cond are final v1 values
%piece together (concatenate) time (t1 and t2) and solution (v1 and v2)
%variables without double counting the seam values
t=[t1; t2(2:end)]; v=[v1; v2(2:end,:)];
%spacetime plot of v variable w/ neurons along y axis, time along x-axis
figure(1); imagesc(v(:,1:N)') ; colorbar
%spacetime plot of r variable w/ neurons along y axis, time along x-axis
figure(2); imagesc(v(:,N+1:end)'); colorbar
%create a movie of the traveling wave
figure(3)
for i=1:length(t)
    plot(v(i,1:N))
    axis([0 N -2.1 2.2])
    pause(.05)
end
function vdot = F_N2(t,v)
global N BC
%
D=.05; a=0.25; b=.001; g=.003; %set the parameters
%
%dV/dt
vdot(1:N^2)=-(v(1:N^2)).*(a-(v(1:N^2))).*(1-(v(1:N^2)))-...
v(N^2+1:end)+D.*secDer2(v(1:N^2),1,BC);
%dR/dt
vdot(N^2+1:2*N^2)=???
%
vdot=vdot';
function V=secDer2(v,dx,BC)
global N
%
F = [0 1 0; 1 -41; 0 1 0]/dx^2;
%determines your boundary conditions
switch BC
case 1 %free bc's
V=conv2(reshape(v,N,N)',F,'same');
V=reshape(V',N*N,1);
case 2 %periodic bc’s
V=conv2_periodic(reshape(v,N,N)',F);
V=reshape(V',N*N,1);
end
function sp = conv2_periodic(s,c)
% 2D convolution for periodic boundary conditions.
% Output of convolution is same size as leading input matrix
[NN,M]=size(s);
[n,m]=size(c); %% both n & m should be odd
%enlarge matrix s in preparation convolution with matrix c via periodic edges
padn = round(n/2) - 1;
padm = round(m/2) - 1;
sp=[zeros(padn,M+(2*padm)); ...
zeros(NN,padm) s zeros(NN,padm); zeros(padn,M+(2*padm))];
%fill in zero padding with the periodic values
sp(1:padn,padm+1:padm+M)=s(NN+1-padn:NN,:);
sp(padn+1+NN:2*padn+NN, padm+1:padm+M)= s(1:padn,:);
sp(padn+1:padn+NN,1:padm)= s(:,M+1-padm:M);
sp(padn+1:padn+NN,padm+M+1:2*padm+M)= s(:,1:padm);
sp(1:padn,1:padm)= s(NN+1-padn:NN,M+1-padm:M);
sp(padn+NN+1:2*padn+NN,1:padm)= s(1:padn,M+1-padm:M);
sp(1:padn,padm+M+1:2*padm+M)=s(NN+1-padn:NN,1:padm);
sp(padn+NN+1:2*padn+NN,padm+M+1:2*padm+M)=s(1:padn,1:padm);
%
%perform 2D convolution
sp = conv2(sp,c,'same');
% reduce matrix back to its original size
sp = sp(padn+1:padn+NN,padm+1:padm+M);
