function isochrons(F,phases,x0)
    % Plot isochrons of a planar dynamical system x’=F(t,x)
    % at points given by the vector phases.
    % x0 is a point on the limit cycle (2x1-vector)
    T = phases(end); % is the period of the cycle
    tau = T/600; % time step of integration
    m = 200; % spatial grid
    k = 5; % the number of skipped cycles
    [t,lc] = ode23s(F,0:tau:T,x0); % forward integration
    dx = (max(lc)-min(lc))’/m; % spatial resolution
    center = (max(lc)+min(lc))’/2; % center of the limit cycle
    iso = [x0-m^0.5*dx, x0+m^0.5*dx]; % isochron’s initial segment
