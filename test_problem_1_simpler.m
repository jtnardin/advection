%test_problem_1.m written 3-30-16 by JTN to run test problem 1 of Thackham
%2009 work. This code is written a bit simpler than test_problem_1.m ,
%which is actually faster and more accurate. Adding matrices (speye +
%A_pos) appears to give some rounding errors that mess up computation, but
%this code should be easier to extend to 2 or 3-D.

n = 50;
dt = 1e-3;

x = linspace(0,1,n);
dx = x(2) - x(1);
t = 0:dt:1;

xn = length(x);
tn = length(t);

x_int = 2:xn-1;

D = 6e-3;
V = -2;
x0 = 0.8;
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

A_pos = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-Dc+-v+v.*sw/2); ...
    (2*Dc+v-v.*se/2-v.*sw/2); (-Dc+v.*se/2)],xn,xn);

A_neg = @(se,sw,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-Dc-v.*sw/2); ...
    (2*Dc+v.*se/2+v.*sw/2-v); (-Dc+v-v.*se/2)],xn,xn);

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
        u(:,i) = (speye(xn) + theta*A_pos(sigma(r_e),sigma(r_w),Vc,x_int,1))\(speye(xn) - (1-theta)*A_pos(sigma(r_e),sigma(r_w),Vc,x_int,1))*u(:,i-1);
        
    elseif V<0
        
        r_e = (u(x_int(1:end-1)+1,i-1) - u(x_int(1:end-1)+2,i-1))./(u(x_int(1:end-1),i-1) - u(x_int(1:end-1)+1,i-1));
        r_e = [r_e;-1];
        r_w = (u(x_int,i-1) - u(x_int+1,i-1))./(u(x_int-1,i-1) - u(x_int,i-1));
        
        %compute interior points
        u(:,i) = (speye(xn) + theta*A_neg(sigma(r_e),sigma(r_w),Vc,x_int,1))\(speye(xn) - (1-theta)*A_neg(sigma(r_e),sigma(r_w),Vc,x_int,1))*u(:,i-1);
        
        %attempt gmres -- only faster for larger systems! (xn = 300, dt = 1e-4)
%         [u(:,i) flag] = gmres((speye(xn) + theta*A_neg(sigma(r_e),sigma(r_w),Vc,x_int,1)),(speye(xn) - (1-theta)*A_neg(sigma(r_e),sigma(r_w),Vc,x_int,1))*u(:,i-1));

        
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

