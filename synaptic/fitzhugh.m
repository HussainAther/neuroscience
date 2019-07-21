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
