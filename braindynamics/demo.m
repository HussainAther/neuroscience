addpath bdtoolkit-2018b/

% construct the model
sys = FHN();

% parameters
%sys.pardef(5).value = 1000;     % Idur = 1000
sys.pardef = bdSetValue(sys.pardef,"Idur",1000);

% time domain
sys.tspan = [0 200];

% solver options
sys.odeoption.RelTol = 1e-5;

% initial conditions
sys.vardef(1).value = 0;    % V=0
sys.vardef(2).value = 0;    % W=0

figure;

for Iamp = linspace(0,5,100)
    % ramp Iamp
    sys.pardef(4).value = Iamp;
    
    % call the solver
    %sol = bdSolve(sys);
    [sys,sol] = bdEvolve(sys,5);

    % extract the results
    t = sol.x;
    V = bdEval(sol,t,1);
    W = bdEval(sol,t,2);

    % plot the results
    %plot(t,V);
    %plot(V,W,".-");
    plot3(Iamp*ones(size(V)),V,W,"k.-");
    title(Iamp);
    xlabel("Iamp");
    ylabel("V");
    zlabel("W");
    xlim([0 5])
    ylim([-3 3]);
    zlim([-1 5]);
    hold on
    
    drawnow;
end

