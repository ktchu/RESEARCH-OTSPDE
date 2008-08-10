%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script computes the solution of the 1d advection equation
%
%    u_t + a u_x = 0
%
% on the domain -\infty < x < \infty with initial conditions
% 
%    u(x) = 0      for x <= -1 
%    u(x) = (x+1)  for -1 < x <= 0
%    u(x) = -(x-1) for 0 < x <= 1
%    u(x) = 0      for x >= 1
%
% using the Lax-Wendroff scheme with a suboptimal time step.  For 
% computational purposes, the domain is truncated to -2 < x < 8 and 
% the boundary conditions are taken to be u(-2) = 0 and u(8) = 0. 
% 
% Kevin T. Chu
% 2007 September
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% PDE parameters
a = 2;

% time integration parameters
t_init  = 0.0;
t_final = 2;

% number of grid cells to use in computation 
% NOTE: number of grid points is N = (N_cell+1)
N_cell = 100;

% construct grid
N = N_cell+1;
dx = 10/N_cell;
x = -2:dx:8;  x = x';   

% set dt
dt = 0.25*dx/a;

% construct central first-derivative operator (with boundary conditions)
G = (diag(ones(N-1,1),1) - diag(ones(N-1,1),-1))/2/dx;
G(1,2) = 0;    % no need to compute gradient at boundary
G(N,N-1) = 0;  % no need to compute gradient at boundary

% construct Laplacian operator (with boundary conditions)
L = (diag(ones(N-1,1),1) - 2*diag(ones(N,1),0) + diag(ones(N-1,1),-1))/dx/dx;
L(1,1) = 0; L(1,2) = 0;   % no need to compute Laplacian at boundary
L(N,N-1) = 0; L(N,N) = 0; % no need to compute Laplacian at boundary

% set initial data
u_init = zeros(size(x));
idx = find(x <= 0 & x > -1);
u_init(idx) = x(idx)+1;
idx = find(x <= 1 & x > 0);
u_init(idx) = -(x(idx)-1);


% first-order upwind (Optimal Time Step)
t = t_init;
u = u_init;
while (t < t_final)

  % compute exact solution
  u_exact = zeros(size(x));
  idx = find((x-a*t) <= 0 & (x-a*t) > -1);
  u_exact(idx) = x(idx)-a*t+1;
  idx = find((x-a*t) <= 1 & (x-a*t) > 0);
  u_exact(idx) = -(x(idx)-a*t-1);

  % plot solution
  figure(1); clf;
  plot(x,u,'b--');
  hold on;
  plot(x,u_exact,'r');
  axis([-2 8 -2 2]);
  title_string = sprintf('t = %f',t);
  title(title_string);
  drawnow

  % impose boundary conditions
  u(1) = 0;
  u(N) = 0;

  % update solution
  u = u - a*dt*(G*u) + a*a*dt^2/2*(L*u);

  % update time
  t = t + dt;

end

% compute exact solution
u_exact = zeros(size(x));
idx = find((x-a*t) <= 0 & (x-a*t) > -1);
u_exact(idx) = x(idx)-a*t+1;
idx = find((x-a*t) <= 1 & (x-a*t) > 0);
u_exact(idx) = -(x(idx)-a*t-1);

% plot final result
figure(1); clf;
plot(x,u,'b--');
hold on;
plot(x,u_exact,'r');
axis([-2 8 -2 2]);
title_string = sprintf('t = %f',t);
title(title_string);

% compute error
err = u-u_exact;
err_L1_upwind1_ots = sum(abs(err))*dx
err_L_inf_upwind1_ots = max((err))

% plot error
figure(2); clf;
plot(x,err,'b-');  


