% FHNode ODE for the canonical FitzHugh-Nagumo model of neural excitability.
%    V' = V - 1/3*V^3 - W + Iapp
%    tau * W' = V + a - b*W
%
% FitzHugh (1961) Impulses and physiological states in theoretical models of nerve membrane. Biophysical J. 1:445?466
% Nagumo, Arimoto and Yoshizawa (1962) An active pulse transmission line simulating nerve axon. Proc. IRE. 50:2061?2070.
function dY = FHNode(~,Y,a,b,tau,Iapp)
