function m = morlet(t, params)
    % Morlet function
    sigma = params(1);
    m = pi^ - 0.25*exp( - i*sigma.*t - 0.=*t.^2);
