%{Response of blowfly photoreceptors to light pulses delivered in the dark (black dots) or from different background 
luminance values (black triangles). The corresponding background luminances are indicated by the arrows on the abscissa. 
The red dots correspond to the steady-state adaptation of the membrane potential to varying background intensities. 
The abscissa is the logarithm of the light intensity (background or background plus pulse) relative to a fixed reference level. 

The ordinate is the membrane potential deflection normalized relative to the maximal deflection, VMax, that could be elicited by a light pulse
Adapted from Laughlin and Hardie's "Common strategies for light adaptation in the peripheral visual systems of fly and dragonfly"
%}

%target membrane potential
v_v = -37.94;

am_v = .1*(25-(v_v+71))./(exp(2.5-(v_v+71)/10)-1);
bm_v = 4*exp(-(v_v+71)/18);
taum_v = 1./(am_v+bm_v);
minf_v = am_v.*taum_v;

ah_v  = 0.07*exp(-(v_v+71)/20);
bh_v = 1./(exp(3-(v_v+71)/10)+1);
tauh_v = 1./(ah_v+bh_v);
hinf_v = ah_v.*tauh_v;

%rates, converted in 1/s = Hz from 1/ms = kHz
am = am_v*1e3;
bm = bm_v*1e3;
ah = ah_v*1e3;
bh = bh_v*1e3;

%these are the mean of the corresponding exponential distributions 
%given by 1/a and 1/b, respectively
m_am = 1/am;
m_bm = 1/bm;

m_ah = 1/ah;
m_bh = 1/bh;

t_sim = 1; %in s, 1000 ms
%simulate m gates. Since average time between transitions is 1.6ms and
%we want 1000 ms of data simulate at least (2000*1.6 = 3200 ms) on average
n_transm = 2000;

%state vector
s_vectm = zeros(3,n_transm);
s_vectm(1:3,1) = rand(3,1) < 0.5;

%transition times
t_vectm = zeros(3,n_transm);

for j = 1:3
    for i = 2: n_transm
        if ( s_vectm(j,i-1) == 0 )
            %closed, duration determined by alpha
            tnext = exprnd(m_am);
            s_vectm(j,i) = 1;
            t_vectm(j,i) = t_vectm(j,i-1) + tnext;
        else
            %open, duration determined by beta
            tnext = exprnd(m_bm);
            s_vectm(j,i) = 0;
            t_vectm(j,i) = t_vectm(j,i-1) + tnext;
        end;
    end;
end;

h_f1 = figure; 
h_a1 = subplot(4,1,1);
h_a2 = subplot(4,1,2);
h_a3 = subplot(4,1,3);
color_vect = ['k' 'k' 'k'];
hdl_vect = [h_a1 h_a2 h_a3];
t_disp = 0.1;
for j = 1:3
    for i = 2:n_transm
        line('Parent',hdl_vect(j),'XData',[t_vectm(j,i-1) t_vectm(j,i)],'YData',[s_vectm(j,i-1) s_vectm(j,i-1)],...
            'Color',color_vect(j));
        line('Parent',hdl_vect(j),'XData',[t_vectm(j,i) t_vectm(j,i)],'YData',[s_vectm(j,i-1) s_vectm(j,i)],...
            'Color',color_vect(j));
    end;
    set(hdl_vect(j),'XLim',[0 t_disp]); %100 ms 
end;

if ( ~isempty(find(t_vectm(:,n_transm) < t_sim)) )
    disp('less than 1000 ms of data for the m-gate, exiting...');
    return;
end;

ylabel(h_a3,'channel state');

%%
%Build vector for m-gates
dt = 0.00005;
t_max = t_sim-dt;
t_vect = 0:dt:t_max;
ind_max = length(t_vect);

s_vm  = zeros(3,ind_max);
for j = 1:3
    i = 2;
    ind_start = 1;
    while  ( t_vectm(j,i) < t_max )
        ind_end = floor(t_vectm(j,i)/dt) + 1;
        s_vm(j,ind_start:ind_end) = s_vectm(j,i-1);
        ind_start = ind_end;
        i = i + 1;
    end;
    s_vm(j,ind_start:ind_max) = s_vectm(j,i-1);
    %figure; plot(t_vect,s_vm(1,:));    
end;

%simulate h gate
n_transh = 100;

s_vecth = zeros(1,n_transh);
s_vecth(1) = rand(1) < 0.5;

t_vecth = zeros(1,n_transh);

for i = 2: n_transh
    if ( s_vecth(i-1) == 0 )
        %closed, duration determined by alpha
        tnext = exprnd(m_ah);
        s_vecth(i) = 1;
        t_vecth(i) = t_vecth(i-1) + tnext;
    else
        %open, duration determined by beta
        tnext = exprnd(m_bh);
        s_vecth(i) = 0;
        t_vecth(i) = t_vecth(i-1) + tnext;
    end;
end;

%figure; 
h_a4 = subplot(4,1,4);
for i = 2:n_transh
    line('Parent',h_a4,'XData',[t_vecth(i-1) t_vecth(i)],'YData',[s_vecth(i-1) s_vecth(i-1)]);
    line('Parent',h_a4,'XData',[t_vecth(i) t_vecth(i)],'YData',[s_vecth(i-1) s_vecth(i)]);
end;
set(h_a4,'XLim',[0 t_disp]); %100 ms 
xlabel(h_a4,'time (s)');

if (t_vecth(n_transh) < t_sim)
    disp('less than 1000 ms of data for the h-gate, exiting...');
    return;
end;

%%

%Build vector for h-gate
s_vh  = zeros(1,ind_max);
i = 2;
ind_start = 1;
while  ( t_vecth(i) < t_max )
    ind_end = floor(t_vecth(i)/dt) + 1;
    s_vh(ind_start:ind_end) = s_vecth(i-1);
    ind_start = ind_end;
    i = i + 1;
end;
s_vh(ind_start:ind_max) = s_vectm(i-1);
%figure; plot(t_vect,s_vm(1,:));    

%compute open probability
%o_v = s_vm(1,:).*s_vm(2,:).*s_vm(3,:).*s_vh(1,:);
%figure;plot(t_vect,o_v);
%line('Parent',handles.axes5,'XData',t_vect,'YData',o_v);
%set(handles.axes5,'XLim',[0 t_sim]);

cum_s = 1 + s_vm(1,:) + s_vm(2,:) + s_vm(3,:) + 4*s_vh(1,:);
%figure; plot(t_vect,cum_s);
h_f2 = figure;
h_a5 = subplot(2,1,1);
line('Parent',h_a5,'XData',t_vect,'YData',cum_s);
set(h_a5,'XLim',[0 2*t_disp]); %200 ms

%%
%carry out single channel simulation using multi-state variable
%build the transition matrix for the entire channel
n_st = 8; %number of states
Qna = zeros(n_st,n_st);

%fill in the non-zero diagonal elements
Qna(1,2) = 3*am;
Qna(1,5) = ah;

Qna(2,1) = bm; 
Qna(2,3) = 2*am;
Qna(2,6) = ah;

Qna(3,2) = 2*bm;
Qna(3,4) = am;
Qna(3,7) = ah;

Qna(4,3) = 3*bm;
Qna(4,8) = ah;

Qna(5,1) = bh;
Qna(5,6) = 3*am;

Qna(6,2) = bh;
Qna(6,5) = bm;
Qna(6,7) = 2*am;

Qna(7,3) = bh;
Qna(7,6) = 2*bm;
Qna(7,8) = am;

Qna(8,4) = bh;
Qna(8,7) = 3*bm;

%sets the diagonal element values as the escape rates from 
%state i
for i = 1:n_st
    Qna(i,i) = sum(Qna(i,:));
end;

%compute the mean escape times
m_escape = 1./diag(Qna);

%cumulative transition probabilities
cQna = zeros(n_st,n_st-1);

%rescale the probabilities to be the conditional
%probabilities given an escape from state i
Qna(1,2:n_st) = Qna(1,2:n_st)/Qna(1,1);
cQna(1,:) = cumsum(Qna(1,2:n_st));

for i = 2:n_st-1
    Qna(i,1:i-1) = Qna(i,1:i-1)/Qna(i,i);
    Qna(i,i+1:n_st) = Qna(i,i+1:n_st)/Qna(i,i);
    cQna(i,:) = cumsum([Qna(i,1:i-1) Qna(i,i+1:n_st)]); 
end;

Qna(n_st,1:n_st-1) = Qna(n_st,1:n_st-1)/Qna(n_st,n_st);
cQna(n_st,:) = cumsum(Qna(n_st,1:n_st-1));

%number of transitions
n_transna = 2000;

s_vectna = zeros(1,n_transna);

%cumulative distribution of the uniform distr. over 8 states
cu_v = [1:8]/8;
%assign equally to each of the possible states
s_vectna(1) = find(rand(1) <= cu_v,1);

t_vectna = zeros(1,n_transna);

for i = 2: n_transna
    c_state = s_vectna(i-1);
    tnext = exprnd(m_escape(c_state));
    t_vectna(i) = t_vectna(i-1) + tnext;
    n_state = find(rand(1)<= cQna(c_state,:),1);
    if ( n_state >= c_state )
        n_state = n_state + 1;
    end;
    s_vectna(i) = n_state;
end;

%figure; handles.axesn = axes;
h_a6 = subplot(2,1,2);
for i = 2:n_transna
    line('Parent',h_a6,'XData',[t_vectna(i-1) t_vectna(i)],'YData',[s_vectna(i-1) s_vectna(i-1)],'Color','r');
    line('Parent',h_a6,'XData',[t_vectna(i) t_vectna(i)],'YData',[s_vectna(i-1) s_vectna(i)],'Color','r');
end;
set(h_a6,'XLim',[0 2*t_disp]); %200 ms 
xlabel(h_a6,'time (s)');

%print(handles.figure1,'-depsc2','nachan.eps');

