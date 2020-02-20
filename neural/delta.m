%  Illustrate the delta rule
r = 0.01;
w1 = 0;
w2 = 0;
eps = 0.01;
for j=1:200,
    cond = ceil(4*rand);
    switch cond
        case 1
            w1 = w1 - r*(w2>1/2+eps);
        case 2
            w2 = w2 - r*(w1>1/2+eps);
        case 3
            w1 = w1 + r*(1-(w1+w2>1/2+eps));
        otherwise
            w2 = w2 + r*(1-(w1+w2>1/2+eps));
    end
    plot(j,w1,"ro")
    hold on
    plot(j,w2,"kx")
end
hold off
legend("w_1","w_2","location","se")
xlabel("iteration","fontsize",14)
ylabel("synaptic weight","fontsize",14)
box off

