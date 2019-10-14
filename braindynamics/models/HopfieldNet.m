% HopfieldNet  Continuous Hopfield network
%   Constructs a continuous Hopfield network with n nodes.
%       tau * V" = -V + W*F(V) + I
%   where V is (nx1) vector on neuron potentials,
%   W is (nxn) matrix of connection weights,
%   I is (nx1) vector of injection currents,
%   F(V)=tanh(b*V) with slope parameter b.
%
% Example:
%   n = 20;                 % number of neurons
%   Wij = 0.5*rand(n);      % random connectivity
%   Wij = Wij + Wij";       % with symmetric connections (Wij=Wji)
%   Wij(1:n+1:end) = 0;     % and non self coupling (zero diagonal)
%   sys = HopfieldNet(Wij); % construct the system struct
%   gui = bdGUI(sys);       % open the Brain Dynamics GUI
function sys = HopfieldNet(Wij)
    % determine the number of nodes from Wij
    n = size(Wij,1);

    % Handle to our ODE function
    sys.odefun = @odefun;
    
    % Our ODE parameters
    sys.pardef = [ struct("name","Wij",  "value",Wij);
                   struct("name","Iapp", "value",rand(n,1));
                   struct("name","b",    "value",1);
                   struct("name","tau",  "value",10) ];
               
    % Our ODE variables
    sys.vardef = struct("name","V", "value",rand(n,1));
    
    % Default time span
    sys.tspan = [0 200];

    % Specify ODE solvers and default options
    sys.odesolver = {@ode45,@ode113,@odeEul};       % ODE solvers
    sys.odeoption = odeset("RelTol",1e-6);          % ODE solver options
 
    % Include the Latex (Equations) panel in the GUI
    sys.panels.bdLatexPanel.title = "Equations"; 
    sys.panels.bdLatexPanel.latex = {"\textbf{HopfieldNet}";
        "";
        "The Continuous Hopfield Network";
        "\qquad $\tau \dot V_i = -V_i + \sum_j W_{ij} \tanh(b\, V_i) + I_{app}$";
        "where";
        "\qquad $V$ is the firing rate of each neuron ($n$ x $1$),";
        "\qquad $K$ is the connectivity matrix ($n$ x $n$),";
        "\qquad $b$ is a slope parameter,";
        "\qquad $I_{app}$ is the applied current ($n$ x $1$),";
        "\qquad $i{=}1 \dots n$.";
        "";
        "Notes";
        ["\qquad 1. This simulation has $n{=}",num2str(n),"$."]};
    
    % Include the Time Portrait panel in the GUI
    sys.panels.bdTimePortrait.title = "Time Portrait";
 
    % Include the Phase Portrait panel in the GUI
    sys.panels.bdPhasePortrait.title = "Phase Portrait";

    % Include the Space-Time panel in the GUI
    sys.panels.bdSpaceTime.title = "Space-Time";

    % Include the Solver panel in the GUI
    sys.panels.bdSolverPanel.title = "Solver";
end

% The ODE function.
function dV = odefun(t,V,Wij,Ii,b,tau)  
    dV = (-V + Wij*tanh(b*V) + Ii)./tau;
end
