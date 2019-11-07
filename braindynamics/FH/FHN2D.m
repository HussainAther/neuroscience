% FHN2D Sheet of FitzHugh-Nagumo neurons
% The spatial equations are
%    V" = V - 1/3*V^3 - W + c^2*(Vxx +Vyy) + Iext
%    tau * W" = V + a - b*W
% where V(x,t) is the membrane voltage at time t for the neuron at
% position x and W(x,t) is the corresponding recovery variable. 
%
% Usage:
%   sys = FHN2D(nr,nc)
% where
%   nr x nc is the number of neurons in the sheet.
%
% Example:
%   nr = 100;                           % number of rows
%   nc = 100;                           % number of columns
%   sys = FHN2D(nr,nc);                 % construct the model
%   gui = bdGUI(sys);                   % run the model in the GUI
%
% Authors
%   Stewart Heitmann (2019)
function sys = FHN2D(nr,nc)
    % Handle to our ODE function
    sys.odefun = @odefun;

    % Default stimulus profile
    rmid = ceil(2*nr/3);
    cmid = ceil(2*nc/3);  
    S = zeros(nr,nc);
    S(rmid-2:rmid+2,cmid-2:cmid+2) = 1;
    
    % Our ODE parameters
    sys.pardef = [
        struct("name","a",     "value",1.0,  "lim",[0 5])
        struct("name","b",     "value",0.5,  "lim",[0 1])
        struct("name","c",     "value",0.5,  "lim",[0 5])
        struct("name","S",     "value",S,    "lim",[0 1])
        struct("name","Iamp",  "value",1,    "lim",[0 5])
        struct("name","Idur",  "value",100,  "lim",[0 50])
        struct("name","tau",   "value",50,   "lim",[1 20])
        struct("name","dx",    "value",1,    "lim",[0.1 10])
        struct("name","dy",    "value",1,    "lim",[0.1 10])
        ];
    
    % Our ODE variables
    sys.vardef = [
        struct("name","V", "value",-1.3*ones(nr,nc), "lim",[-3 3])
        struct("name","W", "value",-0.6*ones(nr,nc), "lim",[-3 3])
        ];
               
    % Default time span
    sys.tspan = [0 10];
    sys.evolve = true;

    % ode options
    sys.odeoption.RelTol = 1e-6;
    sys.odeoption.InitialStep = 0.1;

    % Include the Latex (Equations) panel in the GUI
    sys.panels.bdLatexPanel.title = "Equations"; 
    sys.panels.bdLatexPanel.latex = {
        "\textbf{FitzhughNagumo2D}";
        "";
        num2str([nr nc], "A  %dx%d sheet of FitzHugh-Nagumo neurons");
        "\qquad $\dot V = V - \frac{1}{3}V^3 - W + c^2 (\frac{\partial^2{V}}{\partial x^2} + \frac{\partial^2{V}}{\partial y^2} ) + I$";
        "\qquad $\tau \dot W = V + a - b W$";
        "where";
        "\qquad $V(x,y,t)$ is the membrane voltage at position $x,y$,";
        "\qquad $W(x,y,t)$ is the recovery variable at position $x,y$,";
        "\qquad $I(x,y,t)$ is the injection current at position $x,y$,";
        "\qquad $a,b,c$ are constants,";
        "\qquad $\tau$ is a time constant,";
        "\qquad $dx$ and $dy$ are spatial step sizes.";
        "";
        "The stimulus is defined as";
        "\qquad $I(x,y,t) = I_{amp} \times S(x,y)$";
        "where"
        "\qquad $S(x,y)$ is the spatial profile of the stimulus,";
        "\qquad $Iamp$ is the amplitude of the stimulus,";
        "\qquad $Idur$ is the duration of the stimulus.";
        "";
        "Spatial derivatives are approximated with second-order central differences,";
        "\qquad $\partial^2 V / \partial x^2 \approx \big( V_{i,j-1} - 2V_{i,j} + V_{i,j+1} \big) / dx^2$";
        "\qquad $\partial^2 V / \partial y^2 \approx \big( V_{i-1,j} - 2V_{i,j} + V_{i+1,j} \big) / dy^2$";
        "";
        "\textbf{References}";
        "Fitzhugh (1961) Impulses and physiological states in theoretical models of nerve membrane. Biophysical J. 1:445--466";
        "Nagumo, Arimoto and Yoshizawa (1962) An active pulse transmission line simulating nerve axon. Proc. IRE. 50:2061--2070.";
        "";
...        num2str([n n],"$U$ and $V$ are both %d x %d arrays in this simulation.");
        };
              
    % Other display panels
    sys.panels.bdSpace2D = [];                                   
    sys.panels.bdSolverPanel = [];

              
    % The ODE function. It is implemented as a nested function so that
    % we can inherit variables nr and nc from the enclosing function.
    function dY = odefun(t,Y,a,b,c,S,Iamp,Idur,tau,dx,dy)
        % extract incoming variables
        Y = reshape(Y,[nr nc 2]);   % restore 2D data format
        V = Y(:,:,1);               % 1st plane of Y contains V
        W = Y(:,:,2);               % 2nd plane of Y contains W

        % apply Stimulus
        if (t>=0 && t<Idur)
           Iapp = Iamp*S;
        else
           Iapp = zeros(size(S));
        end

        % Spatial derivatives (with periodic boundary conditions)        
        Vw = circshift(V,-1,2);      % West
        Ve = circshift(V,+1,2);      % East
        Vn = circshift(V,-1,1);      % North
        Vs = circshift(V,+1,1);      % South
        Vxx = (Vw-2*V+Ve) ./dx^2;    % Vxx = (V_{x-1} - 2V_x + V_{x+1} ) / dx^2
        Vyy = (Vs-2*V+Vn) ./dy^2;    % Vyy = (V_{y-1} - 2V_y + V_{y+1} ) / dy^2

        % FitzHugh-Nagumo equations
        dV = V - 1/3*V.^3 - W + c.^2.*(Vxx + Vyy) + Iapp;
        dW = (V + a - b.*W) ./ tau;

        % Return a column vector
        dY = reshape( cat(3,dV,dW), [], 1);
    end

end
