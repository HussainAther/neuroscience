% Continuous wavelet transform (cwt) 
function coefs = simple_cwt(t, x, mother_wavelet, max_wavelet, scales, params)
    % Generates coefs for a continuous wavelet transform
    % t, x are time and data points for time series data
    % mother_wavelet is a function, taking parameters (t, params),
    % where the value of params depends on the specific function used
    % max_wavelet is the maximum range of the wavelet function (beyond which
    % the wavelet is essentially zero)
    % scales is a vector of desired scales
    % params is the parameter for the mother wavelet function
    max_t = max(t);
    dt = t(2)-t(1);
    full_t = - (max_t/2):dt:(max_t/2);
    coefs = zeros(length(scales), length(x));
    points = length(x);
    t_scale = linspace( - max_wavelet, max_wavelet, points);
    dt = (max_wavelet*2)/(points - 1);
    mom_wavelet = feval(mother_wavelet, t_scale, params);
    row = 1;
    for scale = scales
        time_scale = [1 + floor([0:scale*max_wavelet*2]/(scale*dt))];
        wavelet = mom_wavelet(time_scale);
        w = conv(x,wavelet)/sqrt(scale);
        mid_w = floor(length(w)/2);
        mid_x = floor(length(x)/2);
        w = w((( - mid_x:mid_x) + mid_w));
        scale % Print scale to show progress.
        coefs(row,:) = abs(w);
        row = row + 1;
    end 
end
