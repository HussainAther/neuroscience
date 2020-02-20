%{
Upon release, the neurotransmitter molecules enter the synaptic cleft and diffuse across
fairly quickly. On a microscopic level, diffusion is the aggregate effect of many particles
moving randomly. As a first step, examine the motion of a single molecule.
%}

xbounds = [0 10];
ybounds = [0 4];
xdata = [mean(xbounds)];
ydata = [0];
xgrid = 0.01;
ygrid = 0.01;
figure
handle = scatter(xdata, ydata, 'filled');
xlim(xbounds);
ylim(ybounds);

for t = 1:10000
    p = 0.5;
    dx = ((rand > p) - 0.5) * 2;
    dy = ((rand > p) - 0.5) * 2;
    xdata = xdata + dx*xgrid;
    % these two lines assure the molecule stays in x bounds
    xdata(find(xdata < xbounds(1))) = xbounds(1);
    xdata(find(xdata > xbounds(2))) = xbounds(2);
    ydata = ydata + dy*ygrid;
    % these two lines assure the molecule stays within y bounds
    ydata(find(ydata < ybounds(1))) = ybounds(1);
    ydata(find(ydata > ybounds(2))) = ybounds(2);
    set(handle, 'xdata', xdata, 'ydata', ydata);
    drawnow;
end
