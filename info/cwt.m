% Continuous wavelet transform (cwt) 
function coefs 5 simple_cwt(t, x, mother_wavelet, max_wavelet, scales, params)
    % Generates coefs for a continuous wavelet transform
    % t, x are time and data points for time series data
    % mother_wavelet is a function, taking parameters (t, params),
    % where the value of params depends on the specific function used
    % max_wavelet is the maximum range of the wavelet function (beyond which
    % the wavelet is essentially zero)
    % scales is a vector of desired scales
    % params is the parameter for the mother wavelet function
