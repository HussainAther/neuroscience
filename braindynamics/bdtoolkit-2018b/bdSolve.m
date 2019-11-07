%bdSolve  Solve an initial-value problem using the Brain Dynamics Toolbox
%Usage: 
%   sol = bdSolve(sys)
%   sol = bdSolve(sys,tspan)
%   sol = bdSolve(sys,tspan,@solverfun)
%   sol = bdSolve(sys,tspan,@solverfun,solvertype)
%where
%   sys is a system struct describing the dynamical system
%   tspan=[0 100] is the time span of the integration (optional)
%   @solverfun is a function handle to an ode/dde/sde solver (optional)
%   solvertype is a string describing the type of solver (optional).
%
%   The tspan, @solverfun and solvertype arguments are all optional.
%   If tspan is omitted then it defaults to sys.tspan.
%   If @solverfun is omitted then it defaults to the first solver in sys.
%   If @solverfun is supplied but it is not known to the sys struct then
%   you must also supply the solvertype string ("odesolver", "ddesolver"
%   or "sdesolver").
%
%RETURNS
%   sol is the solution structure. It has the same format as that returned
%      by the matlab ode45 solver. Use the bdEval function to extract the
%      results from sol.
%
%EXAMPLE
%   sys = LinearODE;                % Linear system of ODEs
%   tspan = [0 10];                 % integration time domain
%   sol = bdSolve(sys,tspan);       % call the solver
%   tplot = 0:0.1:10;               % time domain of interest
%   Y = bdEval(sol,tplot);          % extract/interpolate the solution
%   plot(tplot,Y);                  % plot the result
%   xlabel("time"); ylabel("y");
function sol = bdSolve(sys,tspan,solverfun,solvertype)
        % check the number of output variables
        if nargout>1
            error("Too many output variables");
        end
        
        % check the number of input variables
        if nargin<1
            error("Not enough input parameters");
        end
   
        % check the validity of the sys struct and fill missing fields with default values
        try
            sys = bd.syscheck(sys);
        catch ME
            throwAsCaller(ME);
        end

        % use defaults for missing input parameters
        switch nargin
            case 1  
                % Caller specified bdSolve(sys).
                
                % Get tspan from sys.tspan 
                tspan = sys.tspan;
                
                % Get solverfun and solvertype from sys 
                if isfield(sys,"odesolver")
                    solverfun = sys.odesolver{1};
                    solvertype = "odesolver";
                end
                if isfield(sys,"ddesolver")
                    solverfun = sys.ddesolver{1};
                    solvertype = "ddesolver";
                end
                if isfield(sys,"sdesolver")
                    solverfun = sys.sdesolver{1};
                    solvertype = "sdesolver";
                end

            case 2
                % Caller specified bdSolve(sys,tspan).
                
                % Get solverfun and solvertype from sys 
                if isfield(sys,"odesolver")
                    solverfun = sys.odesolver{1};
                    solvertype = "odesolver";
                end
                if isfield(sys,"ddesolver")
                    solverfun = sys.ddesolver{1};
                    solvertype = "ddesolver";
                end
                if isfield(sys,"sdesolver")
                    solverfun = sys.sdesolver{1};
                    solvertype = "sdesolver";
                end
                
            case 3
                % Caller specified bdSolve(sys,tspan,solverfun).
                
                % Get solvertype from sys 
                if isfield(sys,"odesolver")
                    solvertype = "odesolver";
                end
                if isfield(sys,"ddesolver")
                    solvertype = "ddesolver";
                end
                if isfield(sys,"sdesolver")
                    solvertype = "sdesolver";
                end
        end
        
        % Call the appropriate solver
        try
            sol = bd.solve(sys,tspan,solverfun,solvertype);
        catch ME
            throwAsCaller(ME);
        end
end
