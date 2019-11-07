doc ode45
Y0 = rand(2,1)
FHNode(0,Y0)
tspan=[0 300]
sol = ode45(@FHNode,tspan,Y0)
sol.x
sol.y
plot(sol.x, sol.y)
plot(sol.x, sol.y(1,:))
plot(sol.y(1,:), sol.y(2,:))
odeset
odeoptions = odeset('RelTol',1e-6)
sol = ode45(@FHNode,tspan,Y0,odeoptions)
plot(sol.y(1,:), sol.y(2,:))
% parameters
a = 1;
b = 0.5;
tau = 10;
Iapp = 1;
whos
sol = ode45(@FHNode,tspan,Y0,odeoptions,a,b,tau,Iapp)
plot(sol.y(1,:), sol.y(2,:),'r')
Iapp
Iapp=0
sol = ode45(@FHNode,tspan,Y0,odeoptions,a,b,tau,Iapp)
plot(sol.y(1,:), sol.y(2,:))
plot(sol.x,sol.y(1,:))
Iapp=0.5;
sol = ode45(@FHNode,tspan,Y0,odeoptions,a,b,tau,Iapp);
plot(sol.x,sol.y(1,:))
Iapp=0.75;
sol = ode45(@FHNode,tspan,Y0,odeoptions,a,b,tau,Iapp);
plot(sol.x,sol.y(1,:))
