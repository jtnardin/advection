%test_problem_2.m written 3-28-16 by JTN to run test problem 2 of Thackham
%2009 work. As of 3-30-16, it is working, though not quite recapitulating
%the results of thackham or pettet. Things to consider in the future
%include:
% - Actually computing V_e and V_w separately -- do not assume constant
% velocity across the CV
% - b(x,t) can probably be computed more accurately. No need to include its
% terms as forcing terms... can include in sparse matrix -- may need flux limiters for these terms, not sure
% - Better incorporate BC in matrix formulation.

%GOING CRAZY!

n = 50;
dt = 1e-4;

x = linspace(0,1,n);
dx = x(2) - x(1);
t = 0:dt:1;

xn = length(x);
tn = length(t);

x_int = 2:xn-1;

theta = 0.5;

%function to access variables for a
a = @(i) xn + i;
%function to access variables for b
b = @(i) 2*xn + i;

%values from Thackham paper
mu_n = 1e-3;
chi = 0.1;
D = 1;
mu_b = 1e-3;
l0 = 100;
l1 = 100;
l2 = 1; %
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

%values from Pettet paper
% mu_n = 1e-3;
% chi = 0.1;
% 
% D = 1;
% mu_b = 1e-3;
% l0 = 100;
% l1 = 100;
% l2 = 10; %
% l3 = 1;
% l4 = 100;
% l5 = 10;
% nhat = 1;
% bhat = 1.5;
% alpha = 2.5;
% l7 = 10;
% xbar = .04;
% bbar = 1.5;
% bstar = 1;
% delta = 1e-2;


%used for computation
Mu_nc = mu_n*dt/dx^2;
Mu_bc = mu_b*dt/dx^2;
Dc = D*dt/dx^2;

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%sparse matrix as a function for computation -- define for positive and
%negative velocity

%positive and negative advection matrices for n
A_np = @(se,sw,v,ind) sparse([ind; ind; ind; 1; xn],[ind-1; ind; ind+1;  1; xn],[(-v*theta+v*theta.*sw/2); ...
    (1+v*theta-v*theta.*se/2-v*theta.*sw/2); (v*theta.*se/2); ones(2,1)],3*xn,3*xn);

A_nm1p = @(se,sw,v,ind) sparse([ind; ind; ind; 1; xn],[ind-1; ind; ind+1;  1; xn],...
    [(Mu_nc+(1-theta)*v-(1-theta)*v.*sw/2); (1-2*Mu_nc-(1-theta)*v+(1-theta)*v.*se/2+(1-theta)*v.*sw/2); ...
    (Mu_nc-(1-theta)*v.*se/2); ones(2,1)],3*xn,3*xn);

A_nn = @(se,sw,v,ind) sparse([ind; ind; ind],[ind-1; ind; ind+1],[(-theta*v.*sw/2); ...
    (1+theta*v.*se/2+theta*v.*sw/2-theta*v); (theta*v-v*theta.*se/2)],3*xn,3*xn);

A_nm1n = @(se,sw,v,ind) sparse([ind; ind; ind],[ind-1; ind; ind+1;],...
    [(Mu_nc+(1-theta)*v.*sw/2); (1-2*Mu_nc+(1-theta)*v-(1-theta)*v.*se/2-(1-theta)*v.*sw/2); ...
    (Mu_nc-(1-theta)*v+(1-theta)*v.*se/2)],3*xn,3*xn);

F_n = @(u,ind) [l1*u(a(ind)).*u(b(ind)) - l2*u(ind) - l0*(u(ind)).^2 ; zeros(2*xn,1)];

% diffusion matrix for a
Da_n = sparse(a(1:xn),a(1:xn),1,3*xn,3*xn);
Da_nm1 = sparse([a(x_int) a(x_int) a(x_int) a(1) a(xn)],[a(x_int-1) a(x_int) a(x_int+1)  a(1) a(xn)],[Dc*ones(1,xn-2) (1-2*Dc)*ones(1,xn-2) Dc*ones(1,xn-2) ones(1,2)],3*xn,3*xn);

F_a = @(u,ind) [zeros(xn,1); l3/2*(1+tanh((bstar-u(b(ind))))) - l4*u(a(ind)) - l5*u(a(ind)).*u(b(ind)) ; zeros(xn,1)];

%diffusion matrix for b -- want to look into more later -- why not include
%all terms in matrix computation, no forcing?

B_n = @(n) sparse([b(x_int) b(x_int) b(x_int) b(1) b(xn)],[b(x_int-1) b(x_int) b(x_int+1) b(1) b(xn)],...
    [-Mu_bc*theta*n(x_int); (1+Mu_bc*theta*(n(x_int)+n(x_int+1))); -Mu_bc*theta*n(x_int+1); ones(2,1)],3*xn,3*xn);

B_nm1 = @(n) sparse([b(x_int) b(x_int) b(x_int) b(1) b(xn)],[b(x_int-1) b(x_int) b(x_int+1)  b(1) b(xn)],...
    [Mu_bc*(1-theta)*n(x_int); (1-Mu_bc*(1-theta)*(n(x_int)+n(x_int+1))); Mu_bc*(1-theta)*n(x_int+1); ones(2,1)],3*xn,3*xn);


%initialize
u = zeros(xn,tn);
u(1:xn,1) = nhat/(xbar^3)*(x-xbar).*(2*x.^2 - xbar*x - xbar^2).*(x<=xbar); %initalize n
u(b(1:xn),1) = (bhat - bbar)/(xbar^3)*(x-xbar).*(2*x.^2 - xbar*x - xbar^2).*(x<=xbar); %initialize b
% a(x,0)=0 already done

%boundary conditions
n_left = @(t) nhat*exp(-alpha*t); %Exponential decay
a_left = @(a1) D/dx/(D/dx+l7*bhat)*a1; %robin BC

for i = 2:tn
    
    %now compute the velocity for n (V = chi*a_x)
    V = chi*central_diff_compute(u(a(1:xn),i-1))/dx;
    
    V = V(2:end-1); %just take interior points for now
    
    %find where the velocity is positive (V_pos) and negative (V_neg)
    V_pos_loc = find(V>=0);
    V_neg_loc = find(V<0);
   
    %find velocity values
    V_pos = V(V_pos_loc);
    V_neg = V(V_neg_loc);
    
    %while not doing boundaries
    V_pos_loc = V_pos_loc + 1;
    V_neg_loc = V_neg_loc + 1;
    
    %compute the appropriate gradient sensors (r_e, r_w) at the locations
    r_e_pos = (u(V_pos_loc,i-1) - u(V_pos_loc-1,i-1))./(u(V_pos_loc+1,i-1) - u(V_pos_loc,i-1));    
    r_w_pos = (u(V_pos_loc(2:end)-1,i-1) - u(V_pos_loc(2:end)-2,i-1))./(u(V_pos_loc(2:end),i-1) - u(V_pos_loc(2:end)-1,i-1));
    r_w_pos = [-1;r_w_pos];
    
    r_e_neg = (u(V_neg_loc+1,i-1) - u(V_neg_loc+2,i-1))./(u(V_neg_loc,i-1) - u(V_neg_loc+1,i-1));
    r_w_neg = (u(V_neg_loc,i-1) - u(V_neg_loc+1,i-1))./(u(V_neg_loc-1,i-1) - u(V_neg_loc,i-1));

    %eliminate NaN values (0/0 -- not steep!)
    r_e_pos(isnan(r_e_pos)) = 1;
    r_w_pos(isnan(r_w_pos)) = 1;
    r_e_neg(isnan(r_e_neg)) = 1;
    r_w_neg(isnan(r_w_neg)) = 1;
    %set inf values to large value
    r_e_pos(isinf(r_e_pos)) = 100;
    r_w_pos(isinf(r_w_pos)) = 100;
    r_e_neg(isinf(r_e_neg)) = 100;
    r_w_neg(isinf(r_w_neg)) = 100;
    
    
    %force term for b
    dndx = central_diff_compute(u((1:xn),i-1))/dx;
    
    F_b = [ zeros(2*xn,1); -mu_n*dndx + chi*u(1:xn,i-1).*central_diff_compute(u(a(1:xn),i-1))/dx];

    A1 = A_np(sigma(r_e_pos),sigma(r_w_pos),V_pos,V_pos_loc) + A_nn(sigma(r_e_neg),sigma(r_w_neg),V_neg,V_neg_loc) + Da_n + B_n(u(1:xn,i-1));
    A2 = A_nm1p(sigma(r_e_pos),sigma(r_w_pos),V_pos,V_pos_loc) + A_nm1n(sigma(r_e_neg),sigma(r_w_neg),V_neg,V_neg_loc) + Da_nm1 + B_nm1(u(1:xn,i-1));
    F = F_n(u(:,i-1),1:xn) + F_a(u(:,i-1),1:xn) + F_b;

    
    u(:,i) = A1\(A2*u(:,i-1)+dt*F);
    

    %boundary conditions
    u(1,i) = n_left(t(i));
    u(b(1),i) = bhat;
    u(a(1),i) = a_left(u(a(2),i));
    
end

for i=1:100:tn
    plot(x,u(1:xn,i),'b')
    hold on
    plot(x,u(a(1:xn),i),'k')
    plot(x,u(b(1:xn),i),'r')
    title(num2str(t(i)));
    axis([0 1 0 2])
    
    legend('n','a','b')
    
    pause(.125)
    hold off
end

