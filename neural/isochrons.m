function isochrons(F,phases,x0)
    % Plot isochrons of a planar dynamical system x’=F(t,x)
    % at points given by the vector phases.
    % x0 is a point on the limit cycle (2x1-vector)
    T = phases(end); % is the period of the cycle
    tau = T/600; % time step of integration
    m = 200; % spatial grid
    k = 5; % the number of skipped cycles
    [t,lc] = ode23s(F,0:tau:T,x0); % forward integration
    dx = (max(lc)-min(lc))'/m; % spatial resolution
    center = (max(lc)+min(lc))’/2; % center of the limit cycle
    iso = [x0-m^0.5*dx, x0+m^0.5*dx]; % isochron’s initial segment
    [t,lc] = ode23s(F,0:tau:T,x0); % forward integration
    dx = (max(lc)-min(lc))'/m; % spatial resolution
    center = (max(lc)+min(lc))’/2; % center of the limit cycle
    iso = [x0-m^0.5*dx, x0+m^0.5*dx]; % isochron's initial segment
    for t=0:-tau:-(k+1)*T % backward integration
        for i=1:size(iso,2)
            iso(:,i)=iso(:,i)-tau*feval(F,t,iso(:,i)); % move one step
        end
        i=1;
        while i<=size(iso,2) % remove infinite solutions
            if any(abs(iso(:,i)-center)>1.5*m*dx) % check boundaries
                iso = [iso(:,1:i-1), iso(:,i+1:end)]; % remove
            else
                i=i+1;
            end;
        end;
        i=1;
        while i<=size(iso,2)-1
            d=sqrt(sum(sum((iso(:,i)-iso(:,i+1))./dx).^2)); % normalized distance
            if d > 2
                iso = [iso(:,1:i), (iso(:,i)+iso(:,i+1))/2 ,iso(:,i+1:end)];
            end;
            if d < 0.5
                iso = [iso(:,1:i), iso(:,i+1:end)];
            else
                i=i+1;
            end;
        end;
        if (mod(-t,T)<=tau/2) & (-t<k*T+tau) % Refresh the screen.
            cla;plot(lc(:,1),lc(:,2),'r'); hold on; % Plot the limit cycle.
        end;
        if min(abs(mod(-t,T)-phases))tau/2 % Plot the isochrons.
            plot(iso(1,:),iso(2,:),'k-'); drawnow;
        end;
    end;
