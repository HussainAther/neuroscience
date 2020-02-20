% Evaluate the size of basins of attraction in hopfield nets.
%
% usage,   hopbasin(100)
%

function hopbasin(N)

close all
linetypek = {'-k' '--k' ':k'};
linetyper = {'-r' '--r' ':r'};

for i=1:2:5,
    
    p = i*N/10;
    
    missfrac1 = zeros(1,N-1);
    missfrac2 = zeros(1,N-1);
    
    for j=1:10,   
        
        S = 2*round(rand(N,p))-1;  % choose random set of p patterns
        missfrac1 = missfrac1 + hopbcheck(S,1);
        missfrac2 = missfrac2 + hopbcheck(S,2);
        
    end
    
    missfrac1 = missfrac1/10;
    missfrac2 = missfrac2/10;
    
    plot(1:N-1,1-missfrac1,linetypek{(i+1)/2})
    hold on
    plot(1:N-1,1-missfrac2,linetyper{(i+1)/2})
    drawnow
    
end

xlabel('overlap with existing memory','fontsize',14)
ylabel('fraction attracted to that memory','fontsize',14)

return

function miss = hopbcheck(S,rule)

W = S*inv(S'*S)*S';
if rule == 2
    W = W - diag(diag(W));
end

[N, p] = size(S);
miss = zeros(1,N-1);

for M=N-1:-1:1,  % overlap size  
    
    s = 10*(1+N-M);
    
    for i=1:s,     % random remainder
        
        pat = 2*round(rand(N,1))-1;
        pat(1:M) = S(1:M,1);
        err = 1;
        iter = 0;
        
        while err > 0 && iter < 20  % flow this pat
            
            y = Hop(W*pat);
            err = norm(y-pat);
            pat = y;
            iter = iter + 1;
            
        end
        
        d = norm(pat-S(:,1));
        if d > 1e-3,
            miss(M) = miss(M) + 1;
        end
        
    end   % i
    
    miss(M) = miss(M)/s;
    
end  % M

return

function y = bin(x)
y = 2*round(x)-1;

function y = Hop(x)
y = (x>0) - (x<=0);
