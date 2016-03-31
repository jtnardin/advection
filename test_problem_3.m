%test_problem_3.m written 3-30-16 by JTN to run test problem 1 of Thackham
%2009 work -- simulate 2 D diffusion-convection problem

%Does not work for n even if x0,y0 = 0.5 (ratiomal?) .... why???

clear all; clc

%initialization, parameters

n = 75;
dt = 1e-3;
x = linspace(0,1,n);
y = linspace(0,1,n);
[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
t = 0:dt:1;
xn = length(x);
yn = length(y);
tn = length(t);

Dx = 1e-2;
Dy = 1e-2;
chix = 8;
chiy = 8;

Dx_c = Dx*dt/dx^2;
Dy_c = Dy*dt/dy^2;

Vx_c = chix*dt/dx;
Vy_c = chiy*dt/dy;

x0 = 0.5;
y0 = 0.5;

theta = 0.5;

%specify boundary nodes
xbd_0 = 1:yn;
xbd_1 = yn*(xn-1)+1:xn*yn;

xbd = union(xbd_0,xbd_1);

ybd_0 = 1:xn:yn*(xn-1)+1;
ybd_1 = yn:yn:xn*yn;

ybd = union(ybd_0,ybd_1);

bd = union(xbd,ybd);

%indices for r_s in computation.
y_ind_1 = 1:yn-2:(xn-2)*(yn-2);
y_ind_nm1 = yn-2:yn-2:(xn-2)*(yn-2);

%specify interior points
xy_int = 1:xn*yn;

xy_int(bd) = [];


exact_soln = @(x,y,t) (1/(4*t+1))*exp(-(x-chix*t-x0).^2/(Dx*(4*t+1))-(y-chiy*t-y0).^2/(Dy*(4*t+1)));

%initial condition
IC = @(x,y) exact_soln(x,y,0);

%boundary conditions

BC_x_0 = @(y,t) exact_soln(0,y,t);
BC_x_1 = @(y,t) exact_soln(1,y,t);
BC_y_0 = @(x,t) exact_soln(x,0,t);
BC_y_1 = @(x,t) exact_soln(x,1,t);

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%initialize
u = zeros(yn*xn,tn);
u0 = IC(X,Y);
u(:,1) = u0(:);


%sparse matrix as a function for computation

A_pos = @(se,sw,D,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-D+-v+v.*sw/2); ...
    (2*D+v-v.*se/2-v.*sw/2); (-D+v.*se/2)],xn*yn,xn*yn);

A_neg = @(se,sw,D,v,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-D-v.*sw/2); ...
    (2*D+v.*se/2+v.*sw/2-v); (-D+v-v.*se/2)],xn*yn,xn*yn);


tic

for i = 2:tn
    
    
    %set BC
    u(xbd_0,i) = BC_x_0(y,t(i));
    u(xbd_1,i) = BC_x_1(y,t(i));
    
    u(ybd_0,i) = BC_y_0(x,t(i));
    u(ybd_1,i) = BC_y_1(x,t(i));
        

    %get Ax matrix
    
    if chix >= 0
        
        r_e = (u(xy_int,i-1) - u(xy_int-yn,i-1))./(u(xy_int+yn,i-1) - u(xy_int,i-1));    
        %for r_w, start at yn-1 to avoid sampling off the grid ... assign
        %value of -1 for first values (would be true for 0 Neumann BC)
        r_w = (u(xy_int(yn-1:end)-yn,i-1) - u(xy_int(yn-1:end)-2*yn,i-1))./(u(xy_int(yn-1:end),i-1) - u(xy_int(yn-1:end)-yn,i-1));
        r_w = [-1*ones(yn-2,1);r_w];
        
        Ax = A_pos(sigma(r_e),sigma(r_w),Dx_c,Vx_c,xy_int,yn);
        
    elseif chix <0
        
        %remove last x strip to avoid sampling off grid ... replace with -1
        r_e = (u(xy_int(1:end-yn+2)+yn,i-1) - u(xy_int(1:end-yn+2)+2*yn,i-1))./(u(xy_int(1:end-yn+2),i-1) - u(xy_int(1:end-yn+2)+yn,i-1));
        r_e = [r_e;-1*ones(yn-2,1)];
        r_w = (u(xy_int,i-1) - u(xy_int+yn,i-1))./(u(xy_int-yn,i-1) - u(xy_int,i-1));
       
        
        Ax = A_neg(sigma(r_e),sigma(r_w),Dx_c,Vx_c,xy_int,yn);
        
    end
    
    
    %get Ay matrix
    
    %get Ax matrix
    
    if chiy >= 0
        
        r_n = (u(xy_int,i-1) - u(xy_int-1,i-1))./(u(xy_int+1,i-1) - u(xy_int,i-1));    
        
        r_s = (u(xy_int-1,i-1) - u(xy_int-2,i-1))./(u(xy_int,i-1) - u(xy_int-1,i-1));
        %for r_s, note that the row corresponding to y=deltay is sampling
        %incorrectly with the second upwind point -- set these sensors = -1
        %(would be true with 0 Neumann BC)
        r_s(y_ind_1) = -1;
        
        Ay = A_pos(sigma(r_n),sigma(r_s),Dy_c,Vy_c,xy_int,1);
        
    elseif chiy < 0
        
        r_n = (u(xy_int+1,i-1) - u(xy_int+2,i-1))./(u(xy_int,i-1) - u(xy_int+1,i-1));
        %for r_n, note that the row corresponding to y=deltay is sampling
        %incorrectly with the second upwind point -- set these sensors = -1
        %(would be true with 0 Neumann BC)
        r_n(y_ind_nm1) = -1;
        r_s = (u(xy_int,i-1) - u(xy_int+1,i-1))./(u(xy_int-1,i-1) - u(xy_int,i-1));
       
        
        Ay = A_neg(sigma(r_n),sigma(r_s),Dy_c,Vy_c,xy_int,1);
        
    end
    
%     (I + theta*(Ap + An))u(:,i) = (I-(1-theta)*(Ap + An))*u(:,i-1)
    
%     u(:,i) = (speye(xn*yn) + theta*(Ax + Ay))\(speye(xn*yn) - (1-theta)*(Ax + Ay))*u(:,i-1);

    [u(:,i),flag] = gmres((speye(xn*yn) + theta*(Ax + Ay)),(speye(xn*yn) - (1-theta)*(Ax + Ay))*u(:,i-1));
      
    
end

toc
% 
for i = 1:10:tn
    
    subplot(1,2,1)
    
    
    contourf(x,y,reshape(u(:,i),yn,xn),'edgecolor','none')
    title('sim')
    view(2)
    
    subplot(1,2,2)
    
    contourf(x,y,exact_soln(X,Y,t(i)),'edgecolor','none')
    title('exact')
    view(2)
    
    pause(.125)
end
