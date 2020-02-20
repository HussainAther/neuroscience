% Evaluate the capacity of a hopfield on a set of hand written numbers from

function hopcapacity

close all
Npat = 150;
N = 400;

[imgs, labels] = readMNIST('train-images.idx3-ubyte','train-labels.idx1-ubyte',Npat,0);
[c1, ia1] = unique(labels);

state = zeros(400,10);
figure(1)

for i=1:10,
    
    subplot(2,5,i)
    imagesc(imgs(:,:,ia1(i)));
    colormap('gray')
    axis equal
    axis off
    
    state(:,i) = bin(reshape(imgs(:,:,ia1(i)),400,1));
    
end

miss1 = hopcheck(state,1);
miss2 = hopcheck(state,2);
miss3 = hopcheck(state,3);

figure(2)
plot((1:10)',miss1,'ks-','linewidth',1.5)
hold on
plot((1:10)',miss2,'ro:','linewidth',1.5)
%plot((1:10)',miss3,'r*--','linewidth',1.5)
hold off
box off
legend('Outer Product','Projection','location','best')
xlabel('number of stored digits','fontsize',14)
ylabel('fraction of misidentified digits','fontsize',14)

% now check convergence of the oja scheme

W = zeros(N);
Wlim = state*inv(state'*state)*state';
I = eye(N);
err = ones(1,1e6);
iter = 0;
M = 10;

while err > 1e-5
    
    i = ceil(M*rand);
    a = state(:,i);
    W = W + (I-W)*(a*a')/N;
    iter = iter + 1;
    err(iter) = norm(W-Wlim,Inf);
    
end

figure(3)
semilogy(1:iter,err(1:iter),'k','linewidth',1.5)
xlabel('iteration, n','fontsize',14)
ylabel('max(|W^n-P(P^TP)^{-1}P^T|)','fontsize',14)
box off
axis tight
return

function miss = hopcheck(state,rule)

miss = zeros(10,1);

for p=2:10,
    
    SP = state(:,1:p);
    switch rule
        case 1
            W = SP*SP';
        case 2
            W = SP*inv(SP'*SP)*SP';
        case 3
            W = SP*inv(SP'*SP)*SP';
            W = W - diag(diag(W));   
    end
    
    for i=1:p
        
        pat = SP(:,i);
        err = 1;
        iter = 0;
        
        while err > 0
            
            y = Hop(W*pat);
            err = norm(y-pat);
            pat = y;
            iter = iter + 1;
            
        end
        
        d = norm(pat-SP(:,i));
        if d>1e-3,
            miss(p)=miss(p)+1;
        end
        
    end   % i
    
    miss(p) = miss(p)/p;
    
end  % p

% figure(3)
% imagesc(W)
%keyboard

return

function y = bin(x)
y = 2*round(x)-1;

function y = Hop(x)
y = (x>=0) - (x<0);
