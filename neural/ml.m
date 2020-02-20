% Find V_r and QA eigenvalues for the Morris-Lecar system (morris lecar)

syms V w

g = struct('Ca',44,'K',80,'L',20);
E = struct('Ca',100,'K',-70,'L',-50);
Cm = 10;

minf = (1 + tanh((V+1)/15))/2;
winf = (1 + tanh(V/30))/2;
tauw = 2/cosh(V/60);

I = g.K*winf*(V-E.K) + g.Ca*minf*(V-E.Ca) + g.L*(V-E.L);
i = inline(char(I));
Vr = fsolve(@(V) i(V), -80)

wr = subs(winf,V,Vr);

F(1) = (winf - w)/tauw;
F(2) = (-g.K*w*(V-E.K) - g.Ca*minf*(V-E.Ca) - g.L*(V-E.L))/Cm;

J = jacobian(F,[w V]);

B = subs(J,{w,V},{wr,Vr})

eig(B)

