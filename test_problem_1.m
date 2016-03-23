%test_problem_1.m written 3-22-16 by JTN to run test problem 1 of Thackham
%2009 work.

%only valid for V >= 0!

n = 100;
dt = 1e-3;

x = linspace(0,1,n);
dx = x(2) - x(1);
t = 0:dt:1;

xn = length(x);
tn = length(t);

x_int = 3:xn-2;


D = 6e-3;
V = 2;
x0 = 0.2;
theta = 0.5;

Dc = D*dt/dx^2;
Vc = V*dt/dx;

%initial condition
IC = @(x) exp(-(x-x0).^2/D);
%left boundary
LB = @(t) 1/sqrt(1+4*t).*exp(-(x0+V*t).^2./(D*(1+4*t)));
%right boundary
RB = @(t) 1/sqrt(1+4*t).*exp(-(1-x0-V*t).^2./(D*(1+4*t)));

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation
% A = @(s_l,s_r) sparse([x_int x_int x_int],[x_int-1 x_int x_int+1],[(Dc + Vc*(-1 + s_r)) ...
%     (1-2*Dc+Vc*(1-s_l/2-s_r/2)) (Dc-Vc*s_l/2)],xn,xn);

A_n = @(se,sw) sparse([x_int x_int x_int  1 2 xn-1 xn],[x_int-1 x_int x_int+1  1 2 xn-1 xn],[(-Vc*theta+Vc*theta*sw/2); ...
    (1+Vc*theta-Vc*theta*se/2-Vc*theta*sw/2); (Vc*theta*se/2); ones(4,1)],xn,xn);

A_nm1 = @(se,sw) sparse([x_int x_int x_int 1 2 xn-1 xn],[x_int-1 x_int x_int+1  1 2 xn-1 xn],...
    [(Dc+(1-theta)*Vc-(1-theta)*Vc*sw/2); (1-2*Dc-(1-theta)*Vc+(1-theta)*Vc*se/2+(1-theta)*Vc*sw/2); ...
    (Dc-(1-theta)*Vc*se/2); ones(4,1)],xn,xn);

%initialize
u = zeros(xn,tn);
u(:,1) = IC(x);

for i = 2:tn
    
    r_e = (u(x_int,i-1) - u(x_int-1,i-1))./(u(x_int+1,i-1) - u(x_int,i-1));    
    r_w = (u(x_int-1,i-1) - u(x_int-2,i-1))./(u(x_int,i-1) - u(x_int-1,i-1));
    
    u(:,i) = A_n(sigma(r_e),sigma(r_w))\A_nm1(sigma(r_e),sigma(r_w))*u(:,i-1);
end

for i = 1:10:tn
    plot(x,u(:,i))
    title(num2str(i));
    axis([0 1 0 1])
    pause(.125)
end

