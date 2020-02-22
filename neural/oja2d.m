%  A two-dimensional example of Oja's rule

%rotation matrix
phi = 30; %degrees
phir = (pi/180)*phi;
rmat = [cos(phir) sin(phir);-sin(phir) cos(phir)];

%make sure we got the right matrix
x = [0; 1];
y = rmat*x;
%figure; plot([0 y(1)],[0 y(2)],'g');
%title('direction of first pca')

%generate 2-d random gaussian variables with 1st pca at 45 deg
%from horizontal
n_iter = 500;
vrndx = normrnd(0,1,1,n_iter);
vrndy = normrnd(0,4,1,n_iter);
vrnd = [vrndx;vrndy];
vrot = rmat*vrnd;
figure(1); plot(vrot(1,:),vrot(2,:),'kx');
hold on;
%plot(4*[0 y(1)],4*[0 y(2)],'r');
%title('random data and direction of pca1');

C = (vrot*vrot')/(n_iter-1);    % The Correlation matrix
[V,D] = eig(C);
d = diag(D);
[md,ind] = max(d);
x = V(:,ind);			% Its principal eigenvector
plot(9*[0 x(1)],9*[0 x(2)],'r');
xlabel('x_{i,1}','fontsize',14)
ylabel('x_{i,2}','fontsize',14)

%intial weight, normalized
w = zeros(2,n_iter);
w(:,1) = rand(2,1);
w(:,1) = w(:,1)/sqrt(w(:,1)'*w(:,1));

figure(2) 

%carry out the iterative oja algorithm
for i = 1:n_iter-1
    a = vrot(:,i);
    b = w(:,i)'*a;
    w(:,i+1) = w(:,i) + (1/(i+200))*b*(a-b*w(:,i));
end
t = 1:n_iter;
plot(t,w(1,:),'k',t,w(2,:),'r')
hold on
xlabel('iteration, i','fontsize',14)
ylabel('synaptic weight','fontsize',14)
legend('w_{i,1}','w_{i,2}','location','best')
box off

%plot final result
%plot(n_iter,w(1),'ko',n_iter,w(2),'ro');

%plot true direction
%plot(n_iter,y(1),'k*',n_iter,y(2),'r*');
%hold off

