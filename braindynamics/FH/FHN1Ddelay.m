% FHN1Ddelay Ring of FitzHugh-Nagumo neurons with transmission delays
%    V" = V - 1/3*V^3 - W + + c^2*Vxx + Iext
%    tau * W" = V + a - b*W
%
% Example:
%   sys = FHN1Ddelay();     % Construct the system struct.
%   gui = bdGUI(sys);       % Open the Brain Dynamics GUI.
%
function sys = FHN1Ddelay(n)    
    % Handle to our DDE function
    sys.ddefun = @FHNdde;
    
    % Our DDE parameters
    sys.pardef = [
        struct("name","a",    "value",1,    "lim",[-5 5])
        struct("name","b",    "value",0.5,  "lim",[-1 1])
        struct("name","c",    "value",0.01, "lim",[0.01 0.05])
        struct("name","S",    "value",zeros(n,1), "lim",[0 1])
        struct("name","tau",  "value",10,   "lim",[1 20])
        struct("name","Iamp", "value",1,    "lim",[0 5])
        struct("name","Idur", "value",1,    "lim",[0 50])
        struct("name","dx",   "value",0.1,  "lim",[0.1 10])
        ];
               
    % Our DDE lag parameters        
    sys.lagdef = [
        struct("name","d", "value",0.01, "lim",[0.01 10])
        ];

    % Our DDE variables        
    sys.vardef = [
        struct("name","V", "value",rand(n,1), "lim",[-3 3])
        struct("name","W", "value",rand(n,1), "lim",[-1 6])        
        ];
        
    % time span
    sys.tspan = [0 100];
    
    % DDE options
    sys.ddesolver = {@dde23,@dde23a};
    sys.ddeoption.RelTol = 1e-6;
    sys.ddeoption.InitialStep = 0.1;
    
    % Equations
    sys.panels.bdLatexPanel.title = "Equations";
    sys.panels.bdLatexPanel.latex = {
        "\textbf{FHN1Ddelay}";
        "";
        "A cable of FitzHugh-Nagumo neurons with transmission delays";
        "\qquad $\dot V = V - \frac{1}{3}V^3 - W + c^2 \frac{\partial^2{V}}{\partial x^2} + I$";
        "\qquad $\tau \dot W = V + a - b W$";
        "where";
        "\qquad $V(x,t)$ is the membrane voltage at position x,";
        "\qquad $W(x,t)$ is the recovery variable at position x,";
        "\qquad $I(x,t)$ is the injection current at position x,";
        "\qquad $a,b,c$ are constants,";
        "\qquad $d$ is the delay between adjacent neurons,";
        "\qquad $\tau$ is a time constant,";
        "\qquad $dx$ is the spatial step between neurons,";
        "";
        "The stimulus is defined as";
        "\qquad $I(x,t) =$ amp $\times$ $S(x)$";
        "where"
        "\qquad $S(x)$ is the spatial profile of the injection current,";
        "\qquad $amp$ modulates the amplitude of the stimulus,";
        "\qquad $onset$ and $dur$ are the onset time and duration of the stimulus.";
        "";
        "The spatial derivative is approximated by the second-order central";
        "central difference,";
        "\qquad $\partial^2 V(t) / \partial x^2 \approx \big( V_{i-1}(t-d) - 2V_{i}(t) + V_{i+1}(t-d) \big) / dx^2$,";
        "with periodic boundary conditions.";
        "";
        "References:";
        "\quad Fitzhugh (1961) Impulses and physiological states in theoretical models of nerve membrane. Biophysical J. 1:445--466";
        "\quad Nagumo, Arimoto and Yoshizawa (1962) An active pulse transmission line simulating nerve axon. Proc. IRE. 50:2061--2070.";
        };
    
    % Display panels
    sys.panels.bdTimePortrait = [];
    sys.panels.bdSpaceTime = [];
    sys.panels.bdSolverPanel = [];

end


% FitzHugh-Nagumo Ordinary Differential Equations
function [dY,Iapp] = FHNdde(t,Y,Z,a,b,c,S,tau,Iamp,Idur,dx)
    % extract incoming variables
    Y = reshape(Y,[],2);
    V = Y(:,1);
    W = Y(:,2);
    
    % extract incoming delayed variables
    Z = reshape(Z,[],2);
    Vz = Z(:,1);
    %Wz = Z(:,2);
    
    % square pulse
    if (t>=0 && t<Idur)
        Iapp = Iamp*S;
    else
        Iapp = 0;
    end
    
    % Spatial Derivative
    Vl = circshift(Vz,-1);
    Vr = circshift(Vz,+1);
    Vxx = (Vl - 2*V + Vr) ./ dx^2;
    
    % system of equations
    dV = V - 1/3*V.^3 - W + c^2*Vxx + Iapp;
    dW = (V + a - b*W) ./ tau;
    
    % return
    dY = [ dV; dW];
end

