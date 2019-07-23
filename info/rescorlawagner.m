%{
Rescorla-Wagner style models describe mechanisms that show how agents
learn to choose the best actions on the basis of rewarding or punishing
feedback. 
}%
nTrials = 100;
lrate = .3;
rewProbs = [.7 .2];
w = .5+zeros(nTrials+1,2);
pPickAct1(ti) = exp(w(ti,1)) / sum(exp(w(ti,:)));
[v1,v2] = deal(-3:.1:3);
for vi=1:length(v2)
    p(:,vi) = exp(v1) ./ (exp(v1)+exp(v2(vi)));
end
imagesc(p)
