%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script computes the solution of the 1d advection equation
%
%    u_t + a u_x = 0
%
% on the domain -\infty < x < \infty with initial conditions
% 
%    u(x) = 1 for x < 0 and u(x) = -1 for x > 0
%
% the Lax-Wendroff scheme with a suboptimal time step.  For computational 
% purposes, the domain is truncated to -1 < x < 1 and the boundary 
% conditions are taken to be u(-1) = 1 and u(1) = -1. 
% 
% Kevin T. Chu
% 2007 September
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% PDE parameters
a = 5;

% time integration parameters
t_init  = 0.0;
t_final = 0.1;

% number of grid cells to use in computation 
% NOTE: number of grid points is N = (N_cell+1)
N_cell = 100;

% construct grid
N = N_cell+1;
dx = 2/N_cell;
x = -1:dx:1;  x = x';   

% set dt
dt = 0.5*dx/a;

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
idx = find(x <= 0);
u_init(idx) = 1;
idx = find(x > 0);
u_init(idx) = -1;


% Lax-Wendroff (suboptimal time step)
t = t_init;
u = u_init;
while (t < t_final)

  % compute exact solution
  u_exact = zeros(size(x));
  idx = find(x <= a*t);
  u_exact(idx) = 1;
  idx = find(x > a*t);
  u_exact(idx) = -1;

  % plot solution
  figure(1); clf;
  plot(x,u,'b--');
  hold on;
  plot(x,u_exact,'r');
  axis([-1 1 -2 2]);
  title_string = sprintf('t = %f',t);
  title(title_string);
  drawnow

  % impose boundary conditions
  u(1) = 1;
  u(N) = -1;

  % update solution
  u = u - a*dt*(G*u) + a*a*dt^2/2*(L*u); 

  % update time
  t = t + dt;

end

% compute exact solution
u_exact = zeros(size(x));
idx = find(x <= a*t);
u_exact(idx) = 1;
idx = find(x > a*t);
u_exact(idx) = -1;

% plot final result
figure(1); clf;
plot(x,u,'b--');
hold on;
plot(x,u_exact,'r');
axis([-1 1 -2 2]);
title_string = sprintf('t = %f',t);
title(title_string);

% compute error
err = u-u_exact;
err_L1_lax_wendroff = sum(abs(err))*dx
err_L_inf_lax_wendroff = max((err))

% plot error
figure(2); clf;
plot(x,err,'b-');  

