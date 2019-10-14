function syscheck(sys)
    %syscheck - Verifies the format of a sys struct.
    %Use syscheck to check the validity of a user-defined sys struct.
    %
    %EXAMPLE
    %  >> sys = LinearODE2D();
    %  >> syscheck(sys);
    %
    %  sys struct format is OK
    %  Calling Y = sys.odefun(t,Y0,a,b,c,d) where
    %  t is size [1  1]
    %  Y0 is size [2  1]
    %  a is size [1  1]
    %  b is size [1  1]
    %  c is size [1  1]
    %  d is size [1  1]
    %  returns Y as size [2 1]
    %  sys.odefun format is OK
    %  ALL TESTS PASSED OK

    % check for deprecated fields that bd.syscheck ordinarily fixes
    % without complaint but we want to want the user to know about them here.
    
    % check for obsolete fields (from version 2017c)
    if isfield(sys,"auxdef") || isfield(sys,"auxfun")
        throwAsCaller(MException("bdtoolkit:syscheck:obsolete","The system structure contains fields relating to auxiliary variables. These fields were deprecated after version 2017c. Remove the ""sys.auxfun"" and ""sys.auxdef"" fields and try again."));               
    end
    if isfield(sys,"self")
        throw(MException("bdtoolkit:syscheck:obsolete","The system structure contains a ""self"" field. This field was deprecated after version 2017c. Remove the ""sys.self"" field and try again."));               
    end

    

    % check that sys is valid and fill missing fields with defaults   
    try
        sys = bd.syscheck(sys);
    catch ME
        throwAsCaller(MException("bdtoolkit:syscheck",ME.message));
    end
    disp("sys struct format is OK");
    
    % test sys.odefun
    if isfield(sys,"odefun")
        if isempty(sys.pardef)
            parnames=",";
        else            
            parnames = sprintf("%s,",sys.pardef.name);
        end
        t = sys.tspan(1);
        Y0 = bdGetValues(sys.vardef);        
        disp(["Calling Y = sys.odefun(t,Y0,",parnames(1:end-1),") where"]);
        disp(["  t is size [", num2str(size(t)), "]"]);
        disp(["  Y0 is size [", num2str(size(Y0)), "]"]);
        for indx = 1:numel(sys.pardef)
            pname = sys.pardef(indx).name;
            psize = size(sys.pardef(indx).value);
            disp(["  " pname " is size [", num2str(psize), "]"]);
        end
        par = {sys.pardef.value};
        Y = sys.odefun(t,Y0,par{:});
        disp(["  ",num2str(size(Y),"Returns Y as size [%d %d]")]);
        if size(Y,2)~=1
            throwAsCaller(MException("bdtoolkit:syscheck","sys.odefun must return Y as a column vector"));
        end
        disp("  sys.odefun format is OK");
    end

    % test sys.odefun with each ODE solver
    if isfield(sys,"odefun")
        if isempty(sys.pardef)
            parnames=",";
        else            
            parnames = sprintf("%s,",sys.pardef.name);
        end
        for solveridx = 1:numel(sys.odesolver)
            solverfun = sys.odesolver{solveridx};
            solverstr = func2str(solverfun);
            disp(["Calling sol = ",solverstr,"(sys.odefun,tspan,Y0,odeoption,",parnames(1:end-1),") returns"]);
            sol = bdSolve(sys,sys.tspan,solverfun);
            solx = num2str(size(sol.x),"[%d %d]");
            soly = num2str(size(sol.y),"[%d %d]");
            disp(["  sol.x is size ", solx]);
            disp(["  sol.y is size ", soly]);
        end
    end
    
    % test sys.ddefun
    if isfield(sys,"ddefun")
        if isempty(sys.pardef)
            parnames=",";
        else
            parnames = sprintf("%s,",sys.pardef.name);
        end
        t = sys.tspan(1);
        Y0 = bdGetValues(sys.vardef);
        lags = bdGetValues(sys.lagdef);
        n = size(Y0,1);
        l = size(lags,1);
        Z = Y0*ones(1,l);
        disp(["Calling Y = sys.ddefun(t,Y0,Z,",parnames(1:end-1),") where"]);
        disp(["  t is size [", num2str(size(t)), "]"]);
        disp(["  Y0 is size [", num2str(size(Y0)), "]"]);
        disp(["  Z is size [", num2str(size(Z)), "]"]);
        for indx = 1:numel(sys.pardef)
            pname = sys.pardef(indx).name;
            psize = size(sys.pardef(indx).value);
            disp(["  " pname " is size [", num2str(psize), "]"]);
        end
        par = {sys.pardef.value};
        Y = sys.ddefun(t,Y0,Z,par{:});
        disp(["  ",num2str(size(Y),"Returns Y as size [%d %d]")]);
        if size(Y,2)~=1
            throwAsCaller(MException("bdtoolkit:syscheck","sys.ddefun must return Y as a column vector"));
        end
        if size(Y,1)~=n
            throwAsCaller(MException("bdtoolkit:syscheck","sys.odefun must return Y as size [%d 1]",n));
        end
        disp("sys.ddefun format is OK");
    end

    % test sys.ddefun and sys.auxfun with each DDE solver
    if isfield(sys,"ddefun")
        if isempty(sys.pardef)
            parnames=",";
        else            
            parnames = sprintf("%s,",sys.pardef.name);
        end
        for solveridx = 1:numel(sys.ddesolver)
            solverfun = sys.ddesolver{solveridx};
            solverstr = func2str(solverfun);
            disp(["Calling sol = ",solverstr,"(sys.ddefun,lags,Y0,tspan,ddeoption,",parnames(1:end-1),") returns"]);
            sol = bdSolve(sys,sys.tspan,solverfun);
            solx = num2str(size(sol.x),"[%d %d]");
            soly = num2str(size(sol.y),"[%d %d]");
            disp(["  sol.x is size ", solx]);
            disp(["  sol.y is size ", soly]);
        end
    end
    
    
    % test sys.sdeF
    if isfield(sys,"sdeF")
        if isempty(sys.pardef)
            parnames=",";
        else            
            parnames = sprintf("%s,",sys.pardef.name);
        end
        t = sys.tspan(1);
        Y0 = bdGetValues(sys.vardef);        
        disp(["Calling Y = sys.sdeF(t,Y0,",parnames(1:end-1),") where"]);
        disp(["  t is size [", num2str(size(t)), "]"]);
        disp(["  Y0 is size [", num2str(size(Y0)), "]"]);
        for indx = 1:numel(sys.pardef)
            pname = sys.pardef(indx).name;
            psize = size(sys.pardef(indx).value);
            disp(["  " pname " is size [", num2str(psize), "]"]);
        end
        par = {sys.pardef.value};
        Y = sys.sdeF(t,Y0,par{:});
        disp(["  ",num2str(size(Y),"Returns Y as size [%d %d]")]);
        if size(Y,2)~=1
            throwAsCaller(MException("bdtoolkit:syscheck","sys.sdeF must return Y as a column vector"));
        end
        disp("  sys.sdeF format is OK");
    end
    
    % test sys.sdeG
    if isfield(sys,"sdeG")
        if isempty(sys.pardef)
            parnames=",";
        else
            parnames = sprintf("%s,",sys.pardef.name);
        end
        t = sys.tspan(1);
        Y0 = bdGetValues(sys.vardef);        
        disp(["Calling G = sys.sdeG(t,Y0,",parnames(1:end-1),") where"]);
        disp(["  t is size [", num2str(size(t)), "]"]);
        disp(["  Y0 is size [", num2str(size(Y0)), "]"]);
        for indx = 1:numel(sys.pardef)
            pname = sys.pardef(indx).name;
            psize = size(sys.pardef(indx).value);
            disp(["  " pname " is size [", num2str(psize), "]"]);
        end
        par = {sys.pardef.value};
        G = sys.sdeG(t,Y0,par{:});
        disp(["  ",num2str(size(G),"Returns G as size [%d %d]")]);      
        if size(G,1)~=size(Y0,1)
            throwAsCaller(MException("bdtoolkit:syscheck","sys.sdeG must return an (nxm) matrix where n=numel(Y0)"));
        end
        if size(G,2)~=sys.sdeoption.NoiseSources
            throwAsCaller(MException("bdtoolkit:syscheck","sys.sdeG must return an (nxm) matrix where m=sys.sdeoption.NoiseSources"));
        end
        disp("  sys.sdeG format is OK");
    end
       
    % test sys.sdeF abd sys.sdeG with each SDE solver
    if isfield(sys,"sdeF") && isfield(sys,"sdeG")
        if isempty(sys.pardef)
            parnames=",";
        else            
            parnames = sprintf("%s,",sys.pardef.name);
        end
        for solveridx = 1:numel(sys.sdesolver)
            solverfun = sys.sdesolver{solveridx};
            solverstr = func2str(solverfun);
            disp(["Calling sol = ",solverstr,"(sys.sdeF,sys.sdeG,tspan,Y0,ddeoption,",parnames(1:end-1),") returns"]);
            sol = bdSolve(sys,sys.tspan,solverfun);
            solx = num2str(size(sol.x),"[%d %d]");
            soly = num2str(size(sol.y),"[%d %d]");
            solW = num2str(size(sol.dW),"[%d %d]");
            disp(["  sol.x is size ", solx]);
            disp(["  sol.y is size ", soly]);
            disp(["  sol.dW is size ", solW]);
        end
    end
  
    disp("ALL TESTS PASSED OK");
end
