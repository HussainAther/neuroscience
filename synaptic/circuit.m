N= 17;

dt = 1; %time step: This should be chosen carefully: remember that the
% in the system is T_fund=RC. If you want to simulate, you should choose
% the time step to be significantly smaller RC, perhaps say, dt=RC/100, see
% below;  but we also have a 10 ms pulse in the system, so let's choose
% dt=0.001
dt=1/1000; % redefinition of dt

V = zeros(10000,N);   %voltage values, go down rows for each time step
%each column is for each V_k
% dV = zeros(size(V)); This is really redundant ( afterwards you can derive
% it from V if you want); just define as a vector as follows:
dV=zeros(1,17);
R = .06;  %I set R=C=1; that is a bit too far from the real time constant 
% in the expeirment, so I reduced C to 0.05 to get closer
C = .1;
t = 0;
Time=linspace(0,10000,17); % I think this is a wrong interpretation of

% linspace. I would define a time axis as after the loop

V(:,1) = 16;  %set V_1 = 16 for each time step
V(:,17) = 0;  %set V_17 = 0 for each time step
dV(:,1) = 0;
dV(:,17) = 0;
%
% the order of operations of the loops should be changed: the inner loop
% (fastest one) should be the spatial part, the outer is the time steps
%
t=dt; % need to initialize
kt=1;
while kt < 10000 % need a discrete index, since the loop refers to
    % matrix elements that depend on the index. To avoid confusion I
    % renamed this to kt
    if t<=10/1000 % need to implement the driving pulse (10ms at 16 Volt)
        V(kt,1)=16;
    else
        V(kt,1)=0;
    end
    for i = 2:16
        % in the following you must make dV(t+1,i) depend on the previous
        % values of V, so the right hand side should have t and NOT t+1; also
        % you cannot use t as a "real" number to index an array, so use k
        % instead
        dV(i) = dt*(V(kt,i-1) - 2*V(kt,i) + V(kt,i+1))/(R*C);%discrete eq
    end
    V(kt+1,:) = V(kt,:) + dV;
    t = t+dt;
    kt=kt+1;
end

Time=(0:dt:t);

% Plot voltage as a function of time.
figure(1);clf
plot(Time, V(:,:))
xlabel('Time (s)')
ylabel('Voltage (V)')
axis([0,.5,0,20])
