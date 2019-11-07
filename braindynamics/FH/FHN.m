% FHN Canonical FitzHugh-Nagumo model of neural excitability
%    V" = V - 1/3*V^3 - W + Iext
%    tau * W" = V + a - b*W
%
% Example:
%   sys = FHN();            % Construct the system struct.
%   gui = bdGUI(sys);       % Open the Brain Dynamics GUI.
function sys = FHN()    
    % Handle to our ODE function
    sys.odefun = @FHNode;
    % Our ODE parameters
    sys.pardef = [
        struct("name","a", "value",1)
        struct("name","b", "value",0.5)
        struct("name","tau", "value",10)
        struct("name","Iapp", "value",1)
        ];
               
    % Our ODE variables        
    sys.vardef = [
        struct("name","V", "value",rand)
        struct("name","W", "value",rand)        
        ];

    % Display panels
    sys.panels = struct();
     
end 
