% Cortical Orientation and Direction Preference Map as
% achieved via a Self Organized Map
% usage:  codpm(G,stims)   e.g.    codpm(128,7e5)

function w = codpm(G,stims)

close all

X = 15; 
Y = 15; 
sigx = 0.5/1;
sigy = 0.5/1;
sigtheta = 0.1;
sigphi = 0.1;

[CX,CY] = meshgrid(1:G,1:G);
rCX = reshape(CX,G*G,1);
rCY = reshape(CY,G*G,1);

w = zeros(G*G,6);
theta = rand(G*G,1)*pi;
phi = rand(G*G,1)*2*pi;
Rtheta0 = 1;     %abs(randn(G*G,1))*sigtheta;
Rphi0 = 1;       %abs(randn(G*G,1))*sigphi;
w(:,1) = Rtheta0.*cos(theta);
w(:,2) = Rtheta0.*sin(theta);
w(:,3) = Rphi0.*cos(phi);
w(:,4) = Rphi0.*sin(phi);

x = zeros(G); 
y = zeros(G); 
for k=1:G
    x(k,:) = k*X/G + (rand(G,1)-1/2)*sigx;
    y(:,k) = k*Y/G + (rand(G,1)-1/2)*sigy;
end
x=x';
y=y';

w(:,5) = reshape(x,G*G,1);
w(:,6) = reshape(y,G*G,1);

wmap(G,w)    % show w

% Train w

eps = 0.02;
sig = 2.5; %50;
sig2 = 2*sig^2;
Rtheta = 1;
Rphi = 1;

v = zeros(1,6);

for t = 1:stims,
    
    theta = rand*pi;                   % generate stimulus
    phi = theta + sign(rand-1/2)*pi/2;
    v(1) = Rtheta*cos(2*theta);
    v(2) = Rtheta*sin(2*theta);
    v(3) = Rphi*cos(phi);
    v(4) = Rphi*sin(phi);
    v(5) = rand*X;
    v(6) = rand*Y;
    
    rv = repmat(v,G*G,1);
    dvw = rv - w;
    err = sum(dvw.^2,2);        % row sums 
    [val,ind] = min(err);

    r2 = (rCX(ind) - rCX).^2 + (rCY(ind) - rCY).^2;
        
    nbd = repmat(eps*exp(-r2/sig2),1,6);   
    w = w + nbd.*dvw;
    
end

wmap(G,w)    % show w

return

function wmap(G,w)

a = reshape(w(:,1),G,G);
b = reshape(w(:,2),G,G);
c = reshape(w(:,3),G,G);
d = reshape(w(:,4),G,G);

figure
Gm = min(G,32);
axis([0 Gm+1 0 Gm+1])
for i=1:Gm,
    xc = i;
    for j=1:Gm,
        yc = Gm-j+1;
        line([xc-a(i,j)/2 xc+a(i,j)/2],[yc-b(i,j)/2 yc+b(i,j)/2],'color','r');
        arrow([xc-c(i,j)/2 yc-d(i,j)/2],[xc+c(i,j)/2 yc+d(i,j)/2],.35,'k');
        hold on
    end
end
axis square
axis off
hold off

x = reshape(w(:,5),G,G);
y = reshape(w(:,6),G,G);
figure
axis off
[C,h] = contour(x,[1:14],'color','r');
hold on
[C,h]=contour(flipud(y),[1:14],'color','b');
axis square
axis off
hold off

return

%
% arrow.m
%
%    arrow(rand(1,2),rand(1,2),.25,'k')
%
function arrow(s,f,a,c)
plot([s(1) f(1)],[s(2) f(2)],c)
hold on
ang = atan2(f(2)-s(2),f(1)-s(1));
r = sqrt((f(1)-s(1))^2+(f(2)-s(2))^2);
z = [s(1)+(1-a)*r*cos(ang) s(2)+(1-a)*r*sin(ang)];
a = a/2;
ze1 = [z(1)+a*r*cos(ang+pi/2) z(2)+a*r*sin(ang+pi/2)];
plot([f(1) ze1(1)],[f(2) ze1(2)],c);
ze1 = [z(1)+a*r*cos(ang-pi/2) z(2)+a*r*sin(ang-pi/2)];
plot([f(1) ze1(1)],[f(2) ze1(2)],c);
hold off
