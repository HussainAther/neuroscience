function N = logistic_map(r,K,tfinal)
    % Simulate the behavior of the Logistic map that is biologically
    % formulated.
    % The Logistic Map used here is
    % N(t+1) = N(t)*exp(r*(1-N/K))

    if(nargin < 3)
        tfinal = 50;
    end
    N = zeros(1,tfinal+1);
    
    N(1) = 1;
    for ii=1:tfinal
        N(ii+1) = N(ii)*exp(r*(1-N(ii)/K));
    end
    
    % Plot the result.
    plot(0:tfinal,N,'o-');
    xlabel('Time');
    ylabel('N');
    title(['Logistic map: r=' num2str(r) ', K=' num2str(K)]);
    
end
