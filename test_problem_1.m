%test_problem_1.m written 3-22-16 by JTN to run test problem 1 of Thackham
%2009 work. This code is harder to understand than test_problem_1_simpler.m
%, but it's actually a but more accurate. Though this code would be hard to
%extend to 2D / 3D.

n = 200;
dt = 1e-3;

x = linspace(0,1,n);
dx = x(2) - x(1);
t = 0:dt:1;

xn = length(x);
tn = length(t);

x_int = 2:xn-1;


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

exact_soln = @(t) 1/sqrt(1+4*t)*exp(-(x-x0-V*t).^2/(D*(1+4*t)));

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation
% A = @(s_l,s_r) sparse([x_int x_int x_int],[x_int-1 x_int x_int+1],[(Dc + Vc*(-1 + s_r)) ...
%     (1-2*Dc+Vc*(1-s_l/2-s_r/2)) (Dc-Vc*s_l/2)],xn,xn);

A_np = @(se,sw) sparse([x_int x_int x_int  1 xn],[x_int-1 x_int x_int+1  1 xn],[(-Vc*theta+Vc*theta*sw/2); ...
    (1+Vc*theta-Vc*theta*se/2-Vc*theta*sw/2); (Vc*theta*se/2); ones(2,1)],xn,xn);

A_nm1p = @(se,sw) sparse([x_int x_int x_int 1 xn],[x_int-1 x_int x_int+1  1 xn],...
    [(Dc+(1-theta)*Vc-(1-theta)*Vc*sw/2); (1-2*Dc-(1-theta)*Vc+(1-theta)*Vc*se/2+(1-theta)*Vc*sw/2); ...
    (Dc-(1-theta)*Vc*se/2); ones(2,1)],xn,xn);

A_nn = @(se,sw) sparse([x_int x_int x_int  1 xn],[x_int-1 x_int x_int+1  1 xn],[(-theta*Vc*sw/2); ...
    (1+theta*Vc*se/2+theta*Vc*sw/2-theta*Vc); (theta*Vc-Vc*theta*se/2); ones(2,1)],xn,xn);

A_nm1n = @(se,sw) sparse([x_int x_int x_int 1 xn],[x_int-1 x_int x_int+1  1 xn],...
    [(Dc+(1-theta)*Vc*sw/2); (1-2*Dc+(1-theta)*Vc-(1-theta)*Vc*se/2-(1-theta)*Vc*sw/2); ...
    (Dc-(1-theta)*Vc+(1-theta)*Vc*se/2); ones(2,1)],xn,xn);



%initialize
u = zeros(xn,tn);
u(:,1) = IC(x);

tic

for i = 2:tn
    
    if V>=0
    
        r_e = (u(x_int,i-1) - u(x_int-1,i-1))./(u(x_int+1,i-1) - u(x_int,i-1));    
        r_w = (u(x_int(2:end)-1,i-1) - u(x_int(2:end)-2,i-1))./(u(x_int(2:end),i-1) - u(x_int(2:end)-1,i-1));
        r_w = [-1;r_w];

        %compute interior points
        u(:,i) = A_np(sigma(r_e),sigma(r_w))\A_nm1p(sigma(r_e),sigma(r_w))*u(:,i-1);
        
    elseif V<0
        
        r_e = (u(x_int(1:end-1)+1,i-1) - u(x_int(1:end-1)+2,i-1))./(u(x_int(1:end-1),i-1) - u(x_int(1:end-1)+1,i-1));
        r_e = [r_e;-1];
        r_w = (u(x_int,i-1) - u(x_int+1,i-1))./(u(x_int-1,i-1) - u(x_int,i-1));
        
        %compute interior points
        u(:,i) = A_nn(sigma(r_e),sigma(r_w))\A_nm1n(sigma(r_e),sigma(r_w))*u(:,i-1);
        
        %attempt gmres -- only faster for larger systems! (xn = 300, dt = 1e-4)
%         [u(:,i) flag] = gmres(A_nn(sigma(r_e),sigma(r_w)),A_nm1n(sigma(r_e),sigma(r_w))*u(:,i-1));

        
    end
    
    %boundary conditions
    u(1,i) = LB(t(i));
    u(end,i) = RB(t(i));
end

toc

for i = 1:10:tn
    plot(x,u(:,i))
    hold on
    plot(x,exact_soln(t(i)),'r')
    title(num2str(i));
    axis([0 1 0 1])
    pause(.125)
    hold off
end

