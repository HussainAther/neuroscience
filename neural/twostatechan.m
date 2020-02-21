% Simulation of a single channel transitioning between closed and open states with rates α = 100 1/s and β = 500 1/s
% simulation of 200 identical channels with 100 of them initially open. The dashed line indicates the average number of 
% open channels at steady state.

%these are the mean of the corresponding exponential distributions 
%given by 1/alpha and 1/beta, respectively
alpha = 100; %in 1/s
m_alpha = 1/alpha;

beta = 500; %in 1/s
m_beta = 1/beta;

n_trans = 30;

s_vect = zeros(1,n_trans);
s_vect(1) = rand(1) < 0.5;

t_vect = zeros(1,n_trans);

for i = 2: n_trans
    if ( s_vect(i-1) == 0 )
        %closed, duration determined by alpha
        tnext = exprnd(m_alpha);
        s_vect(i) = 1;
        t_vect(i) = t_vect(i-1) + tnext;
    else
        %open, duration determined by beta
        tnext = exprnd(m_beta);
        s_vect(i) = 0;
        t_vect(i) = t_vect(i-1) + tnext;
    end;
end;

h_f1 = figure;
h_a1 = subplot(2,1,1);
for i = 2:n_trans
    line('Parent',h_a1,'XData',[t_vect(i-1) t_vect(i)],'YData',[s_vect(i-1) s_vect(i-1)]);
    line('Parent',h_a1,'XData',[t_vect(i) t_vect(i)],'YData',[s_vect(i-1) s_vect(i)]);
end;
set(h_a1,'XLim',[0 100e-3]); %100 ms 

%t_vect(n_trans);

n_chan = 200;
%get approx the same time
n_trans1 = n_trans*n_chan;

%number of open channels
s1_vect = zeros(1,n_trans1);
%s1_vect(1) = round(n_chan*rand(1)); %random number of open channels
s1_vect(1) = 100;

t1_vect = zeros(1,n_trans1);

for i = 2: n_trans1
    %total current transition rate: Nclosed*alpha + Nopen*beta
    nclosed = n_chan-s1_vect(i-1);
    lambda = nclosed*alpha + s1_vect(i-1)*beta;
    
    %time of next transition
    tnext = exprnd(1/lambda);
    t1_vect(i) = t1_vect(i-1) + tnext;
    
    %given a transition, probability of C->O
    lambda1 = nclosed*alpha/lambda;
    
    if ( rand(1) <= lambda1 )
        %C->O transition
        s1_vect(i) = s1_vect(i-1) +1;
    else
        %O->C transition
        s1_vect(i) = s1_vect(i-1) - 1;
    end;
end;

n_max = 1 + 2*(n_trans1-1);
x_v = zeros(1,n_max);
y_v = zeros(1,n_max);
x_v(1) = t1_vect(1);
x_v(2:2:n_max-1) = t1_vect(2:n_trans1);
x_v(3:2:n_max) = t1_vect(2:n_trans1);

y_v(1:2:n_max-2) = s1_vect(1:n_trans1-1);
y_v(2:2:n_max-1) = s1_vect(1:n_trans1-1);
y_v(n_max) = s1_vect(n_trans);

h_a2 = subplot(2,1,2);
line('Parent',h_a2,'XData',x_v,'YData',y_v);

set(h_a2,'XLim',[0 100e-3]); %100 ms 
xlabel('time (s)');
ylabel('number of open channels');

%print(handles.figure1,'-depsc2','twostatechan.eps');

