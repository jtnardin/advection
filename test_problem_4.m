%test_problem_3.m written 4-2-16 by JTN to run test problem 4 of Thackham
%2009 work -- simulate 2 D diffusion-convection problem with
%chemoattractant

%need to fix y=1 boundary.

%should work to simplify process of calculating velocity, location where
%positive and negative velocity occur... external functions?

%Simplify calculation of sensors... definitely should used external
%functions

clear all; clc

%initialization, parameters

n = 45;
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
chix = 0.5;
chiy = 0.5;

l1 = 1e-3;
a0 = 1;
phi0 = 10;

Dx_c = Dx*dt/dx^2;
Dy_c = Dy*dt/dy^2;

Vx_c = chix*dt/dx;
Vy_c = chiy*dt/dy;

x0 = 0.5;
y0 = 0.5;

a = @(x,y) a0*exp(-(x-1/2).^2 - (y-1).^2);

%%%%%%       computing velocities %%%%%%%%%%%%%%%% 
%%%% should this section become its own function? A bit large. %%%%%%%

a_x = diff(a(X,Y),1,1); %x-diff
a_y = diff(a(X,Y),1,2); %y-diff

%compute velocity for interior points
ve = a_x(2:end,2:end-1); %v_w = a_i,j - a_i-1,j. Excluding y boundary points
vw = a_x(1:end-1,2:end-1); %v_e = a_i+1,j - a_i,j

vs =  a_y(2:end-1,1:end-1);%(v_s = a_i,j - a_i,j-1
vn =  a_y(2:end-1,2:end);%v_n = a_i,j+1 - a_i,j

vw = vw(:);
vn = vn(:);
vs = vs(:);
ve = ve(:);

%where are interior velocities positive / negative
Vx_pos_loc = mean([vw ve],2)>= 0;
Vx_neg_loc = mean([vw ve],2) < 0;

Vy_pos_loc = mean([vn vs],2)>= 0;
Vy_neg_loc = mean([vn vs],2) < 0;

%compute velocities at neg, pos points
vw_x_pos = Vx_c*vw(Vx_pos_loc);
ve_x_pos = Vx_c*ve(Vx_pos_loc);

vw_x_neg = Vx_c*vw(Vx_neg_loc);
ve_x_neg = Vx_c*ve(Vx_neg_loc);

vs_y_pos = Vy_c*vs(Vy_pos_loc);
vn_y_pos = Vy_c*vn(Vy_pos_loc);

vs_y_neg = Vy_c*vs(Vy_neg_loc);
vn_y_neg = Vy_c*vn(Vy_neg_loc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%initial condition
IC = @(x,y) phi0*(x==0) + phi0*(x==1) + phi0*(y==0);


%boundary conditions

BC_x_0 = @(y,t) phi0*ones(length(y),1);
BC_x_1 = @(y,t) phi0*ones(length(y),1);
BC_y_0 = @(x,t) phi0*ones(length(x),1);
BC_y_1 = @(x,t) zeros(length(x),1);

%sigma for flux limiters
sigma = @(r) (r+abs(r))./(1+abs(r));

%initialize
u = zeros(yn*xn,tn);
u0 = IC(X,Y);
u0(u0==2*phi0) = phi0;
u(:,1) = u0(:);


%sparse matrix as a function for computation

A_pos = @(se,sw,D,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-D+-vw+vw.*sw/2); ...
    (2*D+ve-ve.*se/2-vw.*sw/2); (-D+ve.*se/2)],xn*yn,xn*yn);

A_neg = @(se,sw,D,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-D-vw.*sw/2); ...
    (2*D+ve.*se/2+vw.*sw/2-vw); (-D+ve-ve.*se/2)],xn*yn,xn*yn);


tic

for i = 2:tn
    
    
    %set BC
    u(xbd_0,i) = BC_x_0(y,t(i));
    u(xbd_1,i) = BC_x_1(y,t(i));
    
    u(ybd_0,i) = BC_y_0(x,t(i));
    u(ybd_1,i) = BC_y_1(x,t(i));
        

    %get Ax matrix by creating positive and negative components
    
    %will need to replace 1st col for Axp and last col for Axn to avoid
    %sampling off grid.
    %locate xy_int points for positive x vel
    xy_intp = xy_int(Vx_pos_loc); 
    %how many points in 1st column
    xy_intp_col1 = numel(find(xy_intp <= 2*yn-1));
    %locate xy_int points for negativ x vel
    xy_intn = xy_int(Vx_neg_loc);
    %how many are in last column
    xy_intn_coln = numel(find(xy_intn >= xn*(yn-2)+2));
        
    %construct Axp (positive matrix)
    r_ep = (u(xy_intp,i-1) - u(xy_intp-yn,i-1))./(u(xy_intp+yn,i-1) - u(xy_intp,i-1));    
    %for r_w, start at yn-1 to avoid sampling off the grid ... assign
    %value of -1 for first values (would be true for 0 Neumann BC)
    r_wp = (u(xy_intp(xy_intp_col1+1:end)-yn,i-1) - u(xy_intp(xy_intp_col1+1:end)-2*yn,i-1))./(u(xy_intp(xy_intp_col1+1:end),i-1) - u(xy_intp(xy_intp_col1+1:end)-yn,i-1));
    r_wp = [-1*ones(xy_intp_col1,1);r_wp];
    
    %replace nan, inf vals
    r_ep(isnan(r_ep)) = 1;
    r_wp(isnan(r_wp)) = 1;
    r_ep(isinf(r_ep)) = 100;
    r_wp(isinf(r_wp)) = 100;
    
    
    %compute A_posx
    Axp = A_pos(sigma(r_ep),sigma(r_wp),Dx_c,ve_x_pos,vw_x_pos,xy_intp,yn);
        

    if  isempty(xy_intn)
        
        %no danger of sampling off grid
        
        r_en = (u(xy_intn+yn,i-1) - u(xy_intn+2*yn,i-1))./(u(xy_intn,i-1) - u(xy_intn+yn,i-1));
        r_wn = (u(xy_intn,i-1) - u(xy_intn+yn,i-1))./(u(xy_intn-yn,i-1) - u(xy_intn,i-1));
 
    else
      
        %remove last x strip to avoid sampling off grid ... replace with -1
        r_en = (u(xy_intn(1:end-xy_intn_coln)+yn,i-1) - u(xy_intn(1:end-xy_intn_coln)+2*yn,i-1))./(u(xy_intn(1:end-xy_intn_coln),i-1) - u(xy_intn(1:end-xy_intn_coln)+yn,i-1));
        r_en = [r_en;-1*ones(xy_intn_coln,1)];
        r_wn = (u(xy_intn,i-1) - u(xy_intn+yn,i-1))./(u(xy_intn-yn,i-1) - u(xy_intn,i-1));


    end
    
    %fix NaN, inf vals
    r_en(isnan(r_en)) = 1;
    r_wn(isnan(r_wn)) = 1;
    r_en(isinf(r_en)) = 100;
    r_wn(isinf(r_wn)) = 100;
    
    
    %compute Axn
    Axn = A_neg(sigma(r_en),sigma(r_wn),Dx_c,ve_x_neg,vw_x_neg,xy_intn,yn);
    
    Ax = Axp + Axn;
            
    %get Ay matrix
    
    %locate xy_int points for positive x vel
    xy_intp = xy_int(Vy_pos_loc); 
    %how many points in 1st row
    xy_intp_row1 = find(mod(xy_intp,yn)==2); %values in row 2
    %locate xy_int points for negativ x vel
    xy_intn = xy_int(Vy_neg_loc);
    %how many are in last column
    xy_intn_rown = find(mod(xy_intn,yn)==yn-1);
       
        
    r_np = (u(xy_intp,i-1) - u(xy_intp-1,i-1))./(u(xy_intp+1,i-1) - u(xy_intp,i-1));    
    r_sp = (u(xy_intp-1,i-1) - u(xy_intp-2,i-1))./(u(xy_intp,i-1) - u(xy_intp-1,i-1));
    %for r_s, note that the row corresponding to y=deltay is sampling
    %incorrectly with the second upwind point -- set these sensors = -1
    %(would be true with 0 Neumann BC)
    r_sp(xy_intp_row1) = -1;

    %correct inf, NaN vals
    r_np(isnan(r_np)) = 1;
    r_sp(isnan(r_sp)) = 1;
    r_np(isinf(r_np)) = 100;
    r_sp(isinf(r_sp)) = 100;
    
    Ayp = A_pos(sigma(r_np),sigma(r_sp),Dy_c,vn_y_pos,vs_y_pos,xy_intp,1);
        
    
        
    r_nn = (u(xy_intn+1,i-1) - u(xy_intn+2,i-1))./(u(xy_intn,i-1) - u(xy_intn+1,i-1));
    %for r_n, note that the row corresponding to y=1-deltay is sampling
    %incorrectly with the second upwind point -- set these sensors = -1
    %(would be true with 0 Neumann BC)
    r_nn(xy_intn_rown) = -1;
    r_sn = (u(xy_intn,i-1) - u(xy_intn+1,i-1))./(u(xy_intn-1,i-1) - u(xy_intn,i-1));

    
    %correct inf, NaN vals
    r_nn(isnan(r_nn)) = 1;
    r_sn(isnan(r_sn)) = 1;
    r_nn(isinf(r_nn)) = 100;
    r_sn(isinf(r_sn)) = 100;
    
    
    Ayn = A_neg(sigma(r_nn),sigma(r_sn),Dy_c,vn_y_neg,vs_y_neg,xy_intn,1);
    
    %remember to incorporate proliferation
    Ay = Ayn + Ayp  - l1*sparse(xy_int,xy_int,1,xn*yn,xn*yn);
        

    [u(:,i),flag] = gmres((speye(xn*yn) + theta*(Ax + Ay)),(speye(xn*yn) - (1-theta)*(Ax + Ay))*u(:,i-1));
      
    
end

toc
% 
for i = 1:10:tn
        
    surf(x,y,reshape(u(:,i),yn,xn),'edgecolor','none')
    view(2)
    
    pause(.125)
    
end
