function plot_cwt(t, coefs, scales)
    imagesc(t, scales, coefs);
    colormap(hot);
    axis xy;
end

% Generate scalogram.
scales = 1:200;
coefs = my_cwt(t, x, @morletfunc, 10, scales, [10]);
plot_cwt(t, coefs, scales);
