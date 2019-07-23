%{
Rescorla-Wagner style models describe mechanisms that show how agents
learn to choose the best actions on the basis of rewarding or punishing
feedback. 
}%
[v1,v2] = deal(-3:.1:3);
for vi=1:length(v2)
    p(:,vi) = exp(v1) ./ (exp(v1)+exp(v2(vi)));
end
imagesc(p)
