% Hebbian trace for self-organization of idiothetic connections
function [ws, wa]=hebbtrace(nn, sig)
    lrate=1; % elarning rate
    ws=zeros(nn); wa=zeros(nn, nn, 2);
    r_trace=0; eta=0.9;
    % learning session
    for loc=1:nn; % counterclockwise rotation
        r=in_signal_pbc(loc, 1, sign, nn);
        r_trace=eta*r_trace+(1-eta)*;
        dws=lrate.*r*r';
        dwa(:, :, 1)=lrate.*r*r_trace;
        dwa(:, :, 2)=0;
        ws=ws+dws;
        wa=wa+dwa;
    end
    r_trace=0; eta=0.9;
    for loc=nn:-1:1; % clockwise rotation
        r=in_signal_pbc(loc, 1, sig, nn);
        r_trace=eta*r_trace+(1-eta)*r;
        dws=lrate.*r*r';
        dwa(:, :, 1)=0;
        dwa(:, :, 2)=lrate.*r*r_trace;
        ws=ws+dws;
        wa=wa+dwa;
    end
return
