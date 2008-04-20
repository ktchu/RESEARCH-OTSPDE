%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script computes the solution of the 1d advection equation
%
%    u_t + a u_x = 0
%
% on the periodic domain -1 < x < 1 with initial conditions
% 
%    u(x) = sin(2*pi*x) + 5*cos(6*pi*x)
%
% using forward Euler time integration with a first-order upwind
% spatial discretization and an optimal time step of dt = dx/a.
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

% visualization parameters
plot_axes = [-1 1 -7 7];

% time integration parameters
t_init  = 0.0;
t_final = 2.0;

% number of grid cells to use in computation
N = 100;  

% construct grid
dx = 2/N;
x = -1:dx:1-dx;  x = x';   % periodic grid point removed
x_plot = [x; x(1)+2];

% set dt 
dt = dx/a;

% construct upwind first-derivative operator (with periodic BCs)
G = (diag(ones(N,1),0) - diag(ones(N-1,1),-1))/dx;
G(1,N) = -1/dx;

% set initial data
u_init = sin(2*pi*x) + 5*cos(6*pi*x);


% first-order upwind (optimal time step)
t = t_init;
u = u_init;
while (t < t_final)

  % plot solution
  figure(1); clf;
  plot(x_plot,[u; u(1)],'bo');
  hold on;
  u_exact = sin(2*pi*(x-a*t)) + 5*cos(6*pi*(x-a*t));
  plot(x_plot,[u_exact; u_exact(1)],'r');
  axis(plot_axes);
  title_string = sprintf('t = %f',t);
  title(title_string);
  drawnow

  % update solution
  u = u - a*dt*(G*u);

  % update time
  t = t + dt;

end

% compute exact solution
u_exact = sin(2*pi*(x-a*t)) + 5*cos(6*pi*(x-a*t));

% plot final result
figure(1); clf;
plot(x_plot,[u; u(1)],'bo');
hold on;
plot(x_plot,[u_exact; u_exact(1)],'r');
axis(plot_axes);
title_string = sprintf('t = %f',t);
title(title_string);

% compute error
err = u-u_exact;
err_L1_upwind1_ots = sum(abs(err))*dx
err_L_inf_upwind1_ots = max((err))

% plot error
figure(2); clf;
plot(x_plot,[err; err(1)],'b-');  

