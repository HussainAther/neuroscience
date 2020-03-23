% Continuous wavelet transform (cwt) 
function coefs 5 simple_cwt(t, x, mother_wavelet, max_wavelet, scales, params)
    % Generates coefs for a continuous wavelet transform
    % t, x are time and data points for time series data
    % mother_wavelet is a function, taking parameters (t, params),
    % where the value of params depends on the specific function used
