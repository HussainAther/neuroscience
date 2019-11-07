% FHN Canonical FitzHugh-Nagumo model of neural excitability
%    V" = V - 1/3*V^3 - W + Iext
%    tau * W" = V + a - b*W
%
% Example:
%   sys = FHN();            % Construct the system struct.
%   gui = bdGUI(sys);       % Open the Brain Dynamics GUI.
%
function sys = FHN()    
    % Handle to our ODE function
    sys.odefun = @FHNode;
    
    % Our ODE parameters
    sys.pardef = [
        struct("name","a",    "value",1,    "lim",[-5 5])
        struct("name","b",    "value",0.5,  "lim",[-1 1])
        struct("name","tau",  "value",10,   "lim",[1 20])
        struct("name","Iapp", "value",1,    "lim",[0 5])
        ];
               
    % Our ODE variables        
    sys.vardef = [
        struct("name","V", "value",rand, "lim",[-3 3])
        struct("name","W", "value",rand, "lim",[-1 6])        
        ];
        
    % time span
    sys.tspan = [0 300];
    
    % ode options
    sys.odeoption.RelTol = 1e-6;
    sys.odeoption.InitialStep = 0.1;
    
    % Equations
    sys.panels.bdLatexPanel.title = "Equations";
    sys.panels.bdLatexPanel.latex = {
        "\textbf{FitzHugh-Nagumo (FHN) Equations}";
        "";
        "The canonical FitzHugh-Nagumo model of neural excitability";
        "\qquad $\dot V = V - \frac{1}{3}V^3 - W + I_{app}$";
        "\qquad $\tau \dot W = V + a - b W$";
        "where";
        "\qquad $V(t)$ is the membrane voltage,";
        "\qquad $W(t)$ is a recovery variable,";
        "\qquad $a,b$ and $\tau$ are constants.";
        "";
        "References:";
        "\quad FitzHugh (1961) Impulses and physiological states in theoretical models of nerve membrane. Biophysical J. 1:445--466";
        "\quad Nagumo, Arimoto and Yoshizawa (1962) An active pulse transmission line simulating nerve axon. Proc. IRE. 50:2061--2070.";
        };
    
    % Display panels
    sys.panels.bdTimePortrait = [];
    sys.panels.bdPhasePortrait = [];
    sys.panels.bdSolverPanel = [];

end

% FitzHugh-Nagumo Ordinary Differential Equations
function dY = FHNode(~,Y,a,b,tau,Iapp)
    % extract incoming variables
    V = Y(1);
    W = Y(2);
    
    % system of equations
    dV = V - 1/3*V^3 - W + Iapp;
    dW = (V + a - b*W) ./ tau;
    
    % return
    dY = [ dV; dW];
end

