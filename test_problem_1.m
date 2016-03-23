%test_problem_1.m written 3-22-16 by JTN to run test problem 1 of Thackham
%2009 work.

%only valid for V >= 0!

n = 25;
dt = 1e-3;

x = linspace(0,1,n);
dx = x(2) - x(1);
t = 0:dt:1;

xn = length(x);
tn = length(t);

x_int = 3:xn-2;


D = 1e-3;
V = 3;
x0 = 0.5;

Dc = D*dt/dx^2;
Vc = D*dt/dx;

%initial condition
IC = @(x) exp(-(x-x0).^2/D);
%left boundary
LB = @(t) 1/sqrt(1+4*t).*exp(-(x0+V*t).^2./(D*(1+4*t)));
%right boundary
RB = @(t) 1/sqrt(1+4*t).*exp(-(1-x0-V*t).^2./(D*(1+4*t)));

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation
A = @(s_l,s_r) sparse([x_int x_int x_int],[x_int-1 x_int x_int+1],[(Dc + Vc*(-1 + s_r)) ...
    (1-2*Dc+Vc*(1-s_l/2-s_r/2)) (Dc-Vc*s_l/2)],xn,xn);

%initialize
u = zeros(xn,tn);
u(:,1) = IC(x);

for i = 2:tn
    
    r_r = (u(x_int,i-1) - u(x_int-1,i-1))./(u(x_int+1,i-1) - u(x_int,i-1));
    
    r_l = (u(x_int-1,i-1) - u(x_int-2,i-1))./(u(x_int,i-1) - u(x_int-1,i-1));
    
    u(:,i) = A(sigma(r_l),sigma(r_r))*u(:,i-1);
end

for i = 1:10:tn
    plot(x,u(:,i))
    title(num2str(i));
    pause(.125)
end

