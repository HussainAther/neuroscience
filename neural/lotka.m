function lotka
    % Lotka-Volterra (lotka volterra) ODEs for showing dynamics of the systems
    % Y1 --[k1]--> 2*Y1            Reproduction of preys
    % Y1 + Y2 --[k2]--> 2*Y2       Prey-Predator interaction
    % Y2 --[k3]-->                 Predators death
    %
    % Y1 is number of prey species, Y2 is  number of predator species
    tspan=[0 100];
    y0=[1;1];
    
    [T,Y]=ode45(@f,tspan,y0);
    figure;
    subplot(121);
    plot(T,Y(:,1),'b*-',T,Y(:,2),'ro-');
    legend('Preys','Predators')
    xlabel('Time')
    ylabel('Y_1, Y_2')
    subplot(122);
    plot(Y(:,1),Y(:,2),'g');
    xlabel('Y_1')
    ylabel('Y_2')
end

function dydt=f(t,y)
    % Kinetic parameters
    k1=1;k2=0.1;k3=0.1;
    % System of ODEs
    dydt = [ k1*y(1)-k2*y(1)*y(2)
    k2*y(1)*y(2)-k3*y(2)];
end
