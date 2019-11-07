% FHN Canonical FitzHugh-Nagumo model of neural excitability
%    V' = V - 1/3*V^3 - W + Iext
%    tau * W' = V + a - b*W
%
% Example:
%   sys = FHN();            % Construct the system struct.
%   gui = bdGUI(sys);       % Open the Brain Dynamics GUI.
