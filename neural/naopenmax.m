%  Compute the membrane potential leading to the 
%  maximum open probability for the sodium channel
%

v = -29:.01:120;
v = v-71;

am = .1*(25-(v+71))./(exp(2.5-(v+71)/10)-1);
bm = 4*exp(-(v+71)/18);
taum = 1./(am+bm);
minf = am.*taum;

ah  = 0.07*exp(-(v+71)/20);
bh = 1./(exp(3-(v+71)/10)+1);
tauh = 1./(ah+bh);
hinf = ah.*tauh;

o_inf = (minf.^3).*hinf;

[p_max ind_m] = max(o_inf);
v_max = v(ind_m);

info_str = sprintf('maximal open probability: %.2g',p_max);
disp(info_str);
info_str = sprintf('corresponding potential: %.4g',v(ind_m));
disp(info_str);

figure; plot(v,o_inf);
hold on;
plot(v_max,p_max,'or');
xlabel('membrane potential');
ylabel('open probability');

