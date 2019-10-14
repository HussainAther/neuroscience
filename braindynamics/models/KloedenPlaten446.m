% KloedenPlaten446 SDE equation (4.46) from Kloeden and Platen (1992)
%   An explicitly solvable Ito SDE from Kloeden and Platen (1992)  
%     dy = -(a + y*b^2)*(1-y^2)*dt + b(1-y^2)*dW
%
% Example:
%   sys = KloedenPlaten446();       % construct the system struct
%   gui = bdGUI(sys);               % open the Brain Dynamics GUI
function sys = KloedenPlaten446()
    % Handles to our SDE functions
    sys.sdeF = @sdeF;               % deterministic part
    sys.sdeG = @sdeG;               % stochastic part
    
    % Our SDE parameters
    sys.pardef = [ struct("name","a", "value",1.0);
                   struct("name","b", "value",0.8) ];

    % Our SDE variables           
    sys.vardef = struct("name","y", "value",0.1);
    
    % Default time span
    sys.tspan = [0 5];
    
    % Specify SDE solvers and default options
    sys.sdesolver = {@sdeEM,@sdeSH};    % Relevant SDE solvers
    sys.sdeoption.InitialStep = 0.005;  % SDE solver step size (optional)
    sys.sdeoption.NoiseSources = 1;     % Number of Wiener noise processes

    % Include the Latex (Equations) panel in the GUI
    sys.panels.bdLatexPanel.title = "Equations"; 
    sys.panels.bdLatexPanel.latex = {"\textbf{KloedenPlaten446}";
        "";
        "Ito Stochastic Differential Equation (4.46) from Kloeden and Platen (1992)";
        "\qquad $dy = -(a + y\,b^2)(1-y^2)\,dt + b(1-y^2)\,dW_t$";
        "where";
        "\qquad $y(t)$ is the dynamic variable,";
        "\qquad $a$ and $b$ are scalar constants.";
        "";
        "It has the explicit solution";
        "\qquad $y = A/B$";
        "where";
        "\qquad $A = (1+y_0) \exp(-2at + 2b W_t) + y_0 - 1$";
        "\qquad $B = (1+y_0) \exp(-2at + 2b W_t) - y_0 + 1$";
        "";
        "";
        "\textbf{Reference}";
        "Kloeden and Platen (1992) Numerical Solution of Stochastic Diffrential Equations"; 
        ""};
        
    % Include the Time Portrait panel in the GUI
    sys.panels.bdTimePortrait = [];

    % Include the Solver panel in the GUI
    sys.panels.bdSolverPanel = [];     
end

% The deterministic part of the equation
function f = sdeF(~,y,a,b)  
    f = -(a + y.*b^2).*(1 - y^2);
end

% The stochastic part of the equation
function G = sdeG(~,y,~,b)  
    G = b.*(1 - y^2);
end
