function sys = SwiftHohenberg1D(n,dx)
    % SwiftHohenberg1D  Discretised Swift-Hohenberg PDE
    %   The Swift-Hohenberg problem
    %        d_t u(x,t) = -(1+d_x^2) u(x,t) - mu*u + nu*u^3 - u^5
    %   discretised in space as a set of coupled ODEs
    %        U' = -(I+Dxx)^2 * U - mu*U + nu*U.^3 - U.^5
    %   where the Laplacian operator Dxx is an nxn matrix
    %        Dxx(i,j) = (U(i+1,j) - 2U(i,j) + U(i+1,j) )/h^2 
    %   and h is the spatial step (dx).
    %
    % Example:
    %   n = 300;                        % number of spatial nodes
    %   dx = 0.25;                      % spatial step size
    %   sys = SwiftHohenberg1D(n,dx);   % construct the system struct
    %   gui = bdGUI(sys);               % open the Brain Dynamics GUI
    % Precompute (part of) the Laplacian, assuming periodic boundary conditions.
    Dxx = sparse( circshift(eye(n),1) -2*eye(n) + circshift(eye(n),-1) );
    % Precompute the identity matrix
    Ix = speye(n);
    % Initial conditions
    x = [1:n]'*dx - (n+1)*dx/2;             % spatial domain, centered on x=0
    U0  = (-tanh(x-8) + tanh(x+8)).*cos(x)+2.0;
    
    % Handle to our ODE function
    sys.odefun = @odefun;
    
    % Our ODE parameters
    sys.pardef = [ struct('name','mu', 'value',1.5);
                   struct('name','nu', 'value',3.0);
                   struct('name','dx', 'value',dx) ];
               
    % Our ODE variables
    sys.vardef = struct('name','U', 'value',U0);
    
    % Default time span
    sys.tspan = [0 20];
         
    % Specify ODE solvers and default options
    %sys.odesolver = {@ode45,@ode23,@ode113,@odeEuler};  % ODE solvers
    sys.odeoption = odeset('RelTol',1e-6);              % ODE solver options

    % Include the Latex (Equations) panel in the GUI
    sys.panels.bdLatexPanel.title = 'Equations'; 
    sys.panels.bdLatexPanel.latex = {'\textbf{SwiftHohenberg1D}';
        '';
        'Spatially discretized Swift-Hohenberg partial differential equation';
        '\qquad $\dot U = -(I + D_{xx})^2 U - \mu U + \nu U^3 - U^5$';
        'where';
        '\qquad $I$ is the Identity matrix,';
        '\qquad $D_{xx}$ is the Laplacian operator,';
        '\qquad $\mu$ and $\nu$ are scalar constants.';
        '';
        'Notes';
        '\qquad 1. $D_{xx,i} = \big( U_{i+1} - 2U_{i} + U_{i+1} \big) / dx^2$.';
        '\qquad 2. $dx$ is the step size of the spatial discretization,';
        '\qquad 3. Boundary conditions are periodic.';
        ['\qquad 4. This simulation has $n{=}',num2str(n),'$.'];
        '';
        'Adapted from Avitabile (2016) Numerical computation of coherent';
        'structures in spatially-extended neural networks. ICMNS 2016.'};
              
    % Include the Time Portrait panel in the GUI
    sys.panels.bdTimePortrait = [];
 
    % Include the Space-Time panel in the GUI
    sys.panels.bdSpaceTime = [];

    % Include the Solver panel in the GUI
    sys.panels.bdSolverPanel = [];                      
              

    % The ODE function; using the precomputed values of Ix and Dxx
    function dU = odefun(t,U,mu,nu,dx)
        % Swift-Hoehenberg equation
        dU = -(Ix + Dxx./dx^2)^2 * U - mu*U + nu*U.^3 - U.^5;
    end

end
  
