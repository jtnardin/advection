%test_problem_2.m written 3-28-16 by JTN to run test problem 2 of Thackham
%2009 work.

n = 150;
dt = 1e-3;

x = linspace(0,1,n);
dx = x(2) - x(1);
t = 0:dt:1;

xn = length(x);
tn = length(t);

x_int = 2:xn-1;

%function to access variables for a
a = @(i) xn + i;
%function to access variables for b
b = @(i) 2*xn + i;

mu_n = 1e-3;
chi = 1.5;
x0 = 0.8;
theta = 1;
D = 1;
mu_b = 1e-3;
l0 = 3;
l1 = 100;
l2 = 1;
l3 = 100;
l4 = 100;
l5 = 10;
nhat = 1;
bhat = 1.5;
alpha = 2.5;
l7 = 10;
xbar = 0.05;
bbar = 0;
bstar = 1;
delta = 0.01;

Mu_nc = mu_n*dt/dx^2;
chic = chi*dt/dx;
Mu_bc = mu_b*dt/dx^2;
Dc = D*dt/dx^2;

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation -- define for positive and
%negative velocity

%for n

A_np = @(se,sw,v,ind) sparse([ind ind ind 1 xn],[ind-1 ind ind+1  1 xn],[(-v*theta+v*theta.*sw/2); ...
    (1+v*theta-v*theta.*se/2-v*theta.*sw/2); (v*theta.*se/2); ones(2,1)],3*xn,3*xn);

A_nm1p = @(se,sw,ind) sparse([ind ind ind 1 xn],[ind-1 ind ind+1  1 xn],...
    [(Mu_nc+(1-theta)*v-(1-theta)*v.*sw/2); (1-2*Mu_nc-(1-theta)*v+(1-theta)*v.*se/2+(1-theta)*v.*sw/2); ...
    (Mu_nc-(1-theta)*v.*se/2); ones(2,1)],3*xn,3*xn);

A_nn = @(se,sw,v,ind) sparse([ind ind ind  1 xn],[ind-1 ind ind+1  1 xn],[(-theta*v.*sw/2); ...
    (1+theta*v.*se/2+theta*v.*sw/2-theta*v); (theta*v-v*theta.*se/2); ones(2,1)],3*xn,3*xn);

A_nm1n = @(se,sw,v,ind) sparse([ind ind ind 1 xn],[ind-1 ind ind+1  1 xn],...
    [(Mu_nc+(1-theta)*v.*sw/2); (1-2*Mu_nc+(1-theta)*v-(1-theta)*v.*se/2-(1-theta)*v.*sw/2); ...
    (Mu_nc-(1-theta)*v+(1-theta)*v.*se/2); ones(2,1)],3*xn,3*xn);

F_n = @(u,ind) l1*u(a(ind)).*u(b(ind)) - l2*u(ind) - l0*(u(ind)).^2;

% for a

Da_n = sparse(a(x_int),a(x_int),1,3*xn,3*xn);
Da_nm1 = sparse([a(x_int) a(x_int) a(x_int)],[a(x_int-1) a(x_int) a(x_int+1)],[Dc*ones(1,xn-2) (1-2*Dc)*ones(1,xn-2) Dc*ones(1,xn-2)],3*xn,3*xn);

F_a = @(u,ind) l3/2*(1+tanh((bstar-u(b(ind))))) - l4*u(a(ind)) - l5*u(a(ind))*u(b(ind));

%for b

B_n = @(n) sparse([b(x_int) b(x_int) b(x_int)],[b(x_int-1) b(x_int) b(x_int+1)],...
    [-Mu_bc*theta*n(x_int) (1+Mu_bc*theta*(n(x_int)+n(x_int+1))) -Mu_bc*theta*n(x_int+1)],3*xn,3*xn);

B_nm1 = @(n) sparse([b(x_int) b(x_int) b(x_int)],[b(x_int-1) b(x_int) b(x_int+1)],...
    [Mu_bc*(1-theta)*n(x_int) (1-Mu_bc*(1-theta)*(n(x_int)+n(x_int+1))) Mu_bc*(1-theta)*n(x_int+1)],3*xn,3*xn);



%initialize
u = zeros(xn,tn);
u(1:xn,1) = nhat/(xbar^3)*(x-xbar).*(2*x.^2 - xbar*x - xbar^2).*(x<=xbar); %initalize n
u(b(1:xn),1) = (bhat - bbar)/(xbar^3)*(x-xbar).*(2*x.^2 - xbar*x - xbar^2).*(x<=xbar); %initialize b
%a(x,0)=0.

for i = 2:tn
    
    %first compute the velocity for n (V = chi*a_x)
    V = chi*central_diff_compute(u(a(1:xn),i-1))/dx;
    
    V = V(2:end-1); %just take interior points for now
    
    %find where the velocity is positive (V_pos) and negative (V_neg)
    V_pos_loc = find(V>=0)+1;
    V_neg_loc = find(V<0)+1;
   
    %find velocity values
    V_pos = V(V_pos_loc);
    V_neg = V(V_neg_loc);
    
    %compute the appropriate gradient sensors (r_e, r_w) at the locations
    r_e_pos = (u(V_pos_loc,i-1) - u(V_pos_loc-1,i-1))./(u(V_pos_loc+1,i-1) - u(V_pos_loc,i-1));    
    r_w_pos = (u(V_pos_loc-1,i-1) - u(V_pos_loc-2,i-1))./(u(V_pos_loc,i-1) - u(x_int(V_pos_loc-1,i-1)));
  
    r_e_neg = (u(V_neg_loc+1,i-1) - u(V_neg_loc+2,i-1))./(u(V_neg_loc,i-1) - u(V_neg_loc+1,i-1));
    r_w_neg = (u(V_neg_loc,i-1) - u(V_neg_loc+1,i-1))./(u(V_neg_loc-1,i-1) - u(V_neg_loc,i-1));

    
    %force term for b
    
    dndx = central_diff_compute(u((1:xn),i-1))/dx;
    dndx = dndx(2:end-1);
    
    F_b = [-mu_n*dndx + V; zeros(2*xn,1)];

%     %compute interior points
%     u(:,i) = (A_np(sigma(r_e_pos),sigma(r_w_pos),V_pos,V_pos) + A_nn(sigma(r_e_neg),sigma(r_w_neg),V_neg,V_neg))...
%         \((A_nm1p(sigma(r_e_pos),sigma(r_w_pos),V_pos,V_pos) + A_nm1n(sigma(r_e_neg),sigma(r_w_neg),V_neg,V_neg))*u(:,i-1) + dt*F_n(u(:,i-1),1:xn));
%         

%     A1*u^n = A2*u^n-1 + dt*F(u) 

    A1 = A_np(sigma(r_e_pos),sigma(r_w_pos),V_pos,V_pos) + A_nn(sigma(r_e_neg),sigma(r_w_neg),V_neg,V_neg) + Da_n + B_n(u(x_int),i-1);
    A2 = A_nm1p(sigma(r_e_pos),sigma(r_w_pos),V_pos,V_pos) + A_nm1n(sigma(r_e_neg),sigma(r_w_neg),V_neg,V_neg) + Da_nm1 + B_nm1(u(x_int,i-1));
    F = F_n(u(:,i-1),1:xn) + F_a(u(:,i-1),1:xn) + F_b;
    
    
%     %boundary conditions
%     u(1,i) = LB(t(i));
%     u(end,i) = RB(t(i));
end

for i = 1:10:tn
    plot(x,u(:,i))
    hold on
    plot(x,exact_soln(t(i)),'r')
    title(num2str(i));
    axis([0 1 0 1])
    pause(.125)
    hold off
end

