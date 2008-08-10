%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script computes the solution of the Korteweg-de Vries (KdV)
% equation
%
%    u_t + u_xxx + 6 u u_x = 0
%
% on the domain -\infty < x < \infty with initial conditions
%
%    u(x) = c/2*( sech(sqrt(c)/2*x) )^2
%
% using forward Euler time integration with a first-order upwind spatial
% discretization and an optimal time step of dt = dx^3/4.  For computational
% purposes, the domain is truncated to x_lo < x < x_hi and the boundary
% conditions are taken to be the appropriate values of the travelling
% wave solution at x = x_lo and x = x_hi.
%
% NOTE:
% - Because we are using higher-order stencils for spatial derivatives
%   we need to impose boundary conditions for a sufficiently wide
%   ghostcell region.
%
% Kevin T. Chu
% 2007 September
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long

% physical parameters
wave_speed = 1;
wave_speed = 10;

% time integration parameters
t_init  = 0.0;
t_final = 1.0;

% number of grid cells to use in computation
% NOTE: number of grid points is N = (N_cell+1)
N_cell = 100;  
N_cell = 200;  
N_cell = 400;  

% set bounds of computational domain
x_lo = -10;
x_hi = 10;

% construct grid
N = N_cell+5;
dx = (x_hi-x_lo)/N_cell;
x = x_lo-2*dx:dx:x_hi+2*dx;  x = x';   % include ghostcells

% set dt 
dt = dx^3/5;  % (suboptimal time step)
dt = dx^3/4;  % (optimal time step)

% visualization parameters
plot_axes = [x_lo x_hi -0.1*wave_speed 0.6*wave_speed];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for main computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct third-order upwind first- and third-derivative operators 
e = ones(N,1);
D1 = 1/dx*spdiags([-1/3*e -1/2*e e -1/6*e], -1:2, N, N); 
D1 = 1/dx*spdiags([1/6*e -e 1/2*e 1/3*e], -2:1, N, N); 


D1 = 1/dx*spdiags([-1/60*e 3/20*e -3/4*e 0*e 3/4*e -3/20*e 1/60*e], ...
                  -3:3, N, N); 
D1 = 1/dx*spdiags([-1/30*e 1/4*e -e 1/3*e 1/2*e -1/20*e], ...
                  -3:2, N, N);   % err = -1/120 u^(5)
D3 = 1/dx^3*spdiags([-7/4*e 25/4*e -17/2*e 11/2*e -7/4*e 1/4*e], -1:4, N, N);
D3 = 1/dx^3*spdiags([-1/4*e -1/4*e 5/2*e -7/2*e 7/4*e -1/4*e], -2:3, N, N); 

% construct fourth-order central second- and fourth-derivative operators 
D2 = 1/dx^2*spdiags([-1/12*e  4/3*e -5/2*e 4/3*e -1/12*e], -2:2, N, N); 
D4 = 1/dx^4*spdiags([-1/6*e 2*e -13/2*e 28/3*e -13/2*e 2*e -1/6*e], ...
                    -3:3, N, N); 

% remove rows for boundary conditions
D1(1,:) = 0; D1(2,:) = 0; D1(3,:) = 0;
D1(N,:) = 0; D1(N-1,:) = 0; D1(N-2,:) = 0;
D2(1,:) = 0; D2(2,:) = 0; D2(3,:) = 0;
D2(N,:) = 0; D2(N-1,:) = 0; D2(N-2,:) = 0;
D3(1,:) = 0; D3(2,:) = 0; D3(3,:) = 0;
D3(N,:) = 0; D3(N-1,:) = 0; D3(N-2,:) = 0;
D4(1,:) = 0; D4(2,:) = 0; D4(3,:) = 0;
D4(N,:) = 0; D4(N-1,:) = 0; D4(N-2,:) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve KdV equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set initial data
u_init = wave_speed/2*(sech(sqrt(wave_speed)/2*x)).^2;

% forward Euler time integration
t = t_init;
time_step = 1;
u = u_init;
while (t < t_final)

  % plot solution
  figure(1); clf;
  plot(x,u,'bo');
  hold on;
  u_exact = wave_speed/2*(sech(sqrt(wave_speed)/2*(x-wave_speed*t))).^2;
  plot(x,u_exact,'r');
  title_string = sprintf('t = %f',t);
  title(title_string);
  cur_axes = axis;
  cur_axes(1) = x_lo;
  cur_axes(2) = x_hi;
  axis(cur_axes);

  % compute error
  err       = u-u_exact;
  err_L1    = norm(err,1)*dx;
  err_L_inf = norm(err,inf);

  % plot error
  figure(2); clf;
  plot(x,err,'b-');  
  title_string = sprintf('t = %f',t);
  title(title_string);
  cur_axes = axis;
  cur_axes(1) = x_lo;
  cur_axes(2) = x_hi;
  axis(cur_axes);
  drawnow

  % display current status of solution
  disp('------------------------------------------------------------------');
  status_str1 = sprintf('Time = %f, Time step = %d', t, time_step);
  status_str2 = sprintf('L inf Error in Solution = %g', err_L_inf);
  status_str3 = sprintf('L1 Error in Solution = %g', err_L1);
  disp(status_str1);
  disp(status_str2);
  disp(status_str3);

  % impose boundary conditions
  u(1) = u_exact(1);
  u(2) = u_exact(2);
  u(3) = u_exact(3);
  u(N) = u_exact(N);
  u(N-1) = u_exact(N-1);
  u(N-2) = u_exact(N-2);

  % update solution
  correction_term = 18*(D2*u).^2 + 30*(D1*u).*(D3*u) ...
                  + 12*u.*(D4*u) + 72*u.*(D1*u).^2 + 36*(u.^2).*(D2*u);
  u = u - dt*(D3*u+6*u.*(D1*u)) + dt^2/2*correction_term;

  % update time and time step
  t = t + dt;
  time_step = time_step + 1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute exact solution
u_exact = wave_speed/2*(sech(sqrt(wave_speed)/2*(x-wave_speed*t))).^2;

% impose boundary conditions
u(1) = u_exact(1);
u(2) = u_exact(2);
u(3) = u_exact(3);
u(N) = u_exact(N);
u(N-1) = u_exact(N-1);
u(N-2) = u_exact(N-2);

% plot final result
figure(1); clf;
plot(x,u,'bo');
hold on;
plot(x,u_exact,'r');
title_string = sprintf('t = %f',t);
title(title_string);
axis(plot_axes);

% compute error
err       = u-u_exact;
err_L1    = norm(err,1)*dx;
err_L_inf = norm(err,inf);

% plot error
figure(2); clf;
plot(x,err,'b-');  
cur_axes = axis;
cur_axes(1) = x_lo;
cur_axes(2) = x_hi;
axis(cur_axes);

% output results statistics
disp(' ');
disp('==================================================================');
disp('Computation Statistics');
err_L_inf_str = sprintf('  L infinity Norm of Error: %g', err_L_inf);
err_L1_str = sprintf('  L1 Norm of Error: %g', err_L1);
disp(err_L_inf_str);
disp(err_L1_str);
disp('==================================================================');


