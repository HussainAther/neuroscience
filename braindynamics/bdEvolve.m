%bdEvolve  Evolve an initial-value problem using the Brain Dynamics Toolbox
%Usage: 
%   [sys,sol] = bdEvolve(sys)
%   [sys,sol] = bdEvolve(sys,rep)
%   [sys,sol] = bdEvolve(sys,rep,tspan)
%   [sys,sol] = bdEvolve(sys,rep,tspan,@solverfun)
%   [sys,sol] = bdEvolve(sys,rep,tspan,@solverfun,solvertype)
%where
%   sys is a system struct describing the dynamical system
%   rep=1 is the number of repeat simulation runs to perform (optional)
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
%   sys is updated with new initial conditions that correspond to the final 
%      conditions of the last simulation run.
%   sol is the solution structure. It has the same format as that returned
%      by the matlab ode45 solver. Use the bdEval function to extract the
%      results from sol.
function [sys,sol] = bdEvolve(sys,rep,tspan,solverfun,solvertype)
        % check the number of output variables
        if nargout>2
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
                % Caller specified bdEvolve(sys).
                
                % Perform one siimulation run only.
                rep = 1;
                
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
                % Caller specified bdEvolve(sys,rep).
                
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
                
            case 3
                % Caller specified bdEvolve(sys,rep,tspan).
                
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
                
            case 4
                % Caller specified bdEvolve(sys,rep,tspan,solverfun).
                
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
        
        try
            % for each repeat simulation
            for r=1:rep
                % Call the appropriate solver
                sol = bd.solve(sys,tspan,solverfun,solvertype);
            
                % Update the initial conditions
                sys.vardef = bdSetValues(sys.vardef,sol.y(:,end));
            end
            
        catch ME
            throwAsCaller(ME);
        end
        
end
