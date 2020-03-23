function plot_cwt(t, coefs, scales)
    imagesc(t, scales, coefs);
    colormap(hot);
    axis xy;
end
