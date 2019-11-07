% FHN Canonical FitzHugh-Nagumo model of neural excitability
%    V" = V - 1/3*V^3 - W + Iext
%    tau * W" = V + a - b*W
%
% Example:
%   sys = FHN();            % Construct the system struct.
%   gui = bdGUI(sys);       % Open the Brain Dynamics GUI.
%
function sys = FHN1D(n)    
    % Handle to our ODE function
    sys.odefun = @FHNode;
    
    % Our ODE parameters
    sys.pardef = [
        struct("name","a",    "value",1,    "lim",[-5 5])
        struct("name","b",    "value",0.5,  "lim",[0 1])
        struct("name","c",    "value",0.01, "lim",[0.01 0.05])
        struct("name","S",    "value",zeros(n,1), "lim",[0 1])
        struct("name","tau",  "value",10,   "lim",[1 20])
        struct("name","Iamp", "value",1,    "lim",[0 5])
        struct("name","Idur", "value",1,    "lim",[0 50])
        struct("name","dx",   "value",0.1,  "lim",[0.1 10])
        ];
               
    % Our ODE variables        
    sys.vardef = [
        struct("name","V", "value",rand(n,1), "lim",[-3 3])
        struct("name","W", "value",rand(n,1), "lim",[-1 6])        
        ];
        
    % time span
    sys.tspan = [0 100];
    
    % ode options
    sys.odeoption.RelTol = 1e-6;
    sys.odeoption.InitialStep = 0.1;
    
    % Equations
    sys.panels.bdLatexPanel.title = "Equations";
    sys.panels.bdLatexPanel.latex = {
        "\textbf{FHN1D}";
        "";
        "A cable of canonical FitzHugh-Nagumo neurons";
        "\qquad $\dot V = V - \frac{1}{3}V^3 - W + c^2 \frac{\partial^2{V}}{\partial x^2} + I$";
        "\qquad $\tau \dot W = V + a - b W$";
        "where";
        "\qquad $V(x,t)$ is the membrane voltage at position x,";
        "\qquad $W(x,t)$ is the recovery variable at position x,";
        "\qquad $I(x,t)$ is the injection current at position x,";
        "\qquad $a,b,c$ are constants,";
        "\qquad $\tau$ is a time constant,";
        "\qquad $dx$ is the spatial step between neurons,";
        "";
        "The stimulus is defined as";
        "\qquad $I(x,t) =$ amp $\times$ $S(x)$";
        "where"
        "\qquad $S(x)$ is the spatial profile of the injection current,";
        "\qquad $I_{amp}$ is the amplitude of the stimulus,";
        "\qquad $I_{dur}$ is the duration of the stimulus.";
        "";
        "The spatial derivative is approximated by the second-order central";
        "central difference,";
        "\qquad $\partial^2 V / \partial x^2 \approx \big( V_{i-1} - 2V_{i} + V_{i+1} \big) / dx^2$,";
        "with periodic boundary conditions.";
        "";
        "References:";
        "\quad Fitzhugh (1961) Impulses and physiological states in theoretical models of nerve membrane. Biophysical J. 1:445--466";
        "\quad Nagumo, Arimoto and Yoshizawa (1962) An active pulse transmission line simulating nerve axon. Proc. IRE. 50:2061--2070.";
        };
    
    % Display panels
    sys.panels.bdTimePortrait = [];
    sys.panels.bdSpaceTime = [];
    sys.panels.bdAuxiliary.auxfun = {@Stimulus};
    sys.panels.bdSolverPanel = [];

end


% FitzHugh-Nagumo Ordinary Differential Equations
function [dY,Iapp] = FHNode(t,Y,a,b,c,S,tau,Iamp,Idur,dx)
    % extract incoming variables
    Y = reshape(Y,[],2);
    V = Y(:,1);
    W = Y(:,2);
    
    % square pulse
    if (t>=0 && t<Idur)
        Iapp = Iamp*S;
    else
        Iapp = 0;
    end
    
    % Spatial Derivative
    Vl = circshift(V,-1);
    Vr = circshift(V,+1);
    Vxx = (Vl - 2*V + Vr) ./ dx^2;
    
    % system of equations
    dV = V - 1/3*V.^3 - W + c^2*Vxx + Iapp;
    dW = (V + a - b*W) ./ tau;
    
    % return
    dY = [ dV; dW];
end


% Auxiliary function for plotting the stimulus profile
function Stimulus(ax,t,sol,a,b,c,S,tau,Iamp,Idur,dx)
    % Reconstruct the stimulus used by FHNode
    n = size(S,1);                  % number of neurons
    tcount = numel(sol.x);          % number of time steps
    Iapp = zeros(n,tcount);
    for idx = 1:tcount
        [~,iapp] = FHNode(sol.x(idx),sol.y(:,idx),a,b,c,S,tau,Iamp,Idur,dx);
        Iapp(:,idx) = iapp;
    end
    
    %plot the stimulus
    pc = pcolor(sol.x,1:n,Iapp);    % pcolor supports variable time steps
    pc.LineStyle = "none";
    pc.Clipping = "off";
    xlim(sol.x([1 end]));
    ylim([1 n]);
    xlabel("time");
    ylabel("space (node)");
    title("Stimulus");
end

