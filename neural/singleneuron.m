%{
Single neuron voltage-gated ion channel.
%}
%These are some default parameter values
I=10;
a=0.02;
b=0.2;
c=-65;
d=8;
%The initial values for v and u
v=-65;
u=b*v;
%Initialize the vector that will contain the membrane potential time series.
v_tot=zeros(1000, 1);
for t=1:1000
    %set v_tot at this time point to the current value of v
    v_tot(t)=v;
    %Reset v and u if v has crossed threshold. See Eq. 3 above.
    if (v>= 30)
        v=c;
        u=u+d;
    end;
    %Use Euler’s method to integrate Eqs. 1 and 2 from above. Here v is
    %calculated in 2 steps in order to keep the time step small (0.5 ms step in the
    %line below).
    v=v+0.5*(0.04*v^2+5*v+140-u+I);
    v=v+0.5*(0.04*v^2+5*vþ140-u+I);
    u=u+a*(b*v-u);
end;
%This line uses the function find to locate the indices of v_tot that hold elements
%with values greater than or equal to 30 and then sets these elements to 30.
%This normalizes to heights of the action potential peaks to 30.
v_tot(find(v_tot >= 30))=30;
%Plot the neuron’s membrane potential.
plot(v_tot);
