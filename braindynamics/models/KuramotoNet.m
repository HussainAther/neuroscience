% KuramotoNet - Network of Kuramoto Phase Oscillators
%   Constructs a Kuramoto network with n nodes.
%       theta_i" = omega_i + SUM_j Kij*sin(theta_i-theta_j)
%   where 
%       theta is an (nx1) vector of oscillator phases (in radians),
%       omega is an (nx1) vector of natural frequencies (cycles/sec)
%       Kij is an (nxn) matrix of connection weights,
%
% Example:
%   n = 20;                    % number of oscillators
%   Kij = ones(n);             % coupling matrix
%   sys = KuramotoNet(Kij);    % construct the system struct
%   gui = bdGUI(sys);          % open the Brain Dynamics GUI
function sys = KuramotoNet(Kij)
    % Determine the number of nodes from the size of the coupling matrix
    n = size(Kij,1);
    
    % Handle to our ODE function
    sys.odefun = @odefun;
    
    % ODE parameters
    sys.pardef = [ struct("name","Kij",   "value",Kij);
                   struct("name","k",     "value",1);               
                   struct("name","omega", "value",randn(n,1)) ];
               
    % ODE state variables
    sys.vardef = struct("name","theta", "value",2*pi*rand(n,1), "lim",[-pi pi]);
    
    % Time span
    sys.tspan = [0 100];

    % Relevant ODE solvers
    sys.odesolver = {@ode45,@ode23,@ode113,@odeEul};
    
    % ODE solver options
    sys.odeoption.RelTol = 1e-6;
    sys.odeoption.MaxStep = 0.1;                
                    
    % Latex (Equations) panel
    sys.panels.bdLatexPanel.title = "Equations"; 
    sys.panels.bdLatexPanel.latex = {"\textbf{KuramotoNet}";
        "";
        "A generalised network of Kuramoto Oscillators";
        "\qquad $\dot \theta_i = \omega_i + \frac{k}{n} \sum_j K_{ij} \sin(\theta_i - \theta_j)$";
        "where";
        "\qquad $\theta_i$ is the phase of the $i^{th}$ oscillator (radians),";
        "\qquad $\omega_i$ is its natural oscillation frequency (cycles/sec),";
        "\qquad $K$ is the network connectivity matrix ($n$ x $n$),";
        "\qquad $k$ is a scaling constant,";
        "\qquad $i,j=1 \dots n,$";
        ["\qquad $n{=}",num2str(n),".$"];
        "";
        "The Kuramoto order parameter ($R$) is a metric of phase synchronisation.";
        "\qquad $R = \frac{1}{n} \| \sum_i \exp(\mathbf{i} \theta_i) \|$";
        "It corresponds to the radius of the centroid of the phases, as shown in";
        "the Auxiliary panel.";
        "";
        "References";
        "\qquad Kuramoto (1984) Chemical oscillations, waves and turbulence.";
        "\qquad Strogatz (2000) From Kuramoto to Crawford.";
        "\qquad Breakspear et al (2010) Generative models of cortical oscillations.";
        "\qquad Chapter 6.2 of the Handbook for the Brain Dynamics Toolbox (Version 2018b)."};
    
    % Time Portrait panel
    sys.panels.bdTimePortrait.title = "Time Portrait";
    sys.panels.bdTimePortrait.mod = true;
 
    % Phase Portrait panel
    sys.panels.bdPhasePortrait.title = "Phase Portrait";
    sys.panels.bdPhasePortrait.mod = true;

    % Auxiliary panel
    sys.panels.bdAuxiliary.title = "Auxiliary";
    sys.panels.bdAuxiliary.auxfun = {@centroid1,@centroid2,@KuramotoR};
    
    % Solver panel
    sys.panels.bdSolverPanel.title = "Solver";                
end

% Kuramoto ODE function where
% theta is a (1xn) vector of oscillator phases,
% Kij is either a scalar or an (nxn) matrix of connection weights,
% k is a scalar,
% omega is either a scalar or (1xn) vector of oscillator frequencies.
function dtheta = odefun(t,theta,Kij,k,omega)
    n = numel(theta);
    theta_i = theta * ones(1,n);                        % (nxn) matrix with same theta values in each row
    theta_j = ones(n,1) * theta";                       % (nxn) matrix with same theta values in each col
    theta_ij = theta_i - theta_j;                       % (nxn) matrix of all possible (theta_i - theta_j) combinations
    dtheta = omega + k/n.*sum(Kij.*sin(theta_ij),1)";   % Kuramoto Equation in vector form.
end

% Auxiliary function that plots the centroid of the oscillators
function centroid1(ax,t,sol,Kij,k,omega)
    % Get the phases of the oscillators at time t
    theta = bdEval(sol,t);
    
    % Project the phases into the complex plane.
    ztheta = exp(1i.*theta);
    
    % Plot the centroid.
    centroidplot(ztheta);
    text(-1,-1,num2str(t,"time = %g"));
    title("centroid of oscillators"); 
end

% Auxiliary function that plots the centroid of the oscillators
% in a rotating frame where the first oscillator is pinned at
% zero phase.
function centroid2(ax,t,sol,Kij,k,omega)
    % Get the phases of the oscillators at time t
    theta = bdEval(sol,t);
    
    % Project the phases into the complex plane and
    % rotate the frame to pin the first oscillator at zero phase.
    ztheta = exp(1i.*theta) .* exp(-1i*theta(1));
    
    % Plot the centroid.
    centroidplot(ztheta);
    text(-1,-1,num2str(t,"time = %g"));
    title("centroid (rotating frame)"); 
end

% Utility function for plotting the centroid
function centroidplot(ztheta)
    % compute the phase centroid
    centroid = mean(ztheta);
    
    % plot the unit circle
    plot(exp(1i.*linspace(-pi,pi,100)), "color",[0.75 0.75 0.75]);
    
    % plot the oscillator phases on the unit circle
    plot(ztheta,"o","color","k");
    
    % plot the centroid (yellow paddle)
    plot([0 centroid], "color", "k");
    plot(centroid,"o","color","k", "Marker","o", "MarkerFaceColor","y", "MarkerSize",10);
    
    % axis limits etc
    axis equal;
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
end

% Auxiliary function for plotting the Kuramoto order parameter (R).
function KuramotoR(ax,t,sol,Kij,k,omega)
    % Project the phases into the complex plane.
    ztheta = exp(1i.*sol.y);

    % compute the running phase centroid
    centroid = mean(ztheta);

    % plot the amplitide of the centroid versus time.
    %axis normal;
    plot(sol.x,abs(centroid),"color","k","linewidth",1.5);
    
    % axis limits etc
    t0 = sol.x(1);
    t1 = sol.x(end);
    xlim([t0 t1]);
    ylim([-0.1 1.1]);
    xlabel("time");
    ylabel("R = abs(centroid)");
    title("Kuramoto Order Parameter (R)");
end
