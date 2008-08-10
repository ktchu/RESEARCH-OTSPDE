%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for the 2d advection
% equation.
%  
% Kevin T. Chu
% 2008 February
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact
set(0,'DefaultAxesFontSize',18,'DefaultAxesFontName','Helvetica')
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultTextFontSize',18,'DefaultTextFontName','Helvetica')

% set print format
use_color_figures = 0;
print_suffix = 'eps';
if use_color_figures
  print_format = 'epsc';
else
  print_format = 'eps';
end
fig_dir = 'figures';
if ~exist(fig_dir, 'dir')
  mkdir(fig_dir);
end

% set flag for loading data from saved files (instead of recomputing solution)
use_saved_data = 0;
data_dir = 'data-adv_eqn_2d_no_src';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% set simulation parameters
debug_on  = 0;
timing_on = 1;

% physical parameters
A_x = -1.0;  % flow speed in x-direction
A_y = -2.0;  % flow speed in y-direction

% time integration parameters
% NOTE: t_init = 0.0
t_final = 3.0;

% start clock for timing plot generation time
t_start = cputime;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100;

% set dx and dt
dx = 20/N;
dt_FE = 0.5*dx/( abs(A_x) + abs(A_y) );  % suboptimal time step

% solve advection equation using forward Euler with OTS
disp('---------------------------');
disp('  Forward Euler OTS')
disp('---------------------------');
[u_FE_OTS, u_exact, dy_FE_OTS, X_FE_OTS, Y_FE_OTS, timing_data_FE_OTS] = ...
   solveAdvectionEqnForwardEulerOTS2d(A_x, A_y, ...
                                      dx, ...
                                      t_final, ...
                                      debug_on, timing_on);

% solve diffusion equation using forward Euler
disp('---------------------------');
disp('  Forward Euler')
disp('---------------------------');
[u_FE, u_exact, X_FE, Y_FE, timing_data_FE] = ...
   solveAdvectionEqnForwardEuler2d(A_x, A_y, ...
                                   dx, dt_FE, ...
                                   t_final, ...
                                   debug_on, timing_on);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_plot = 101;
x_lo = -10; x_hi = 10;
dx_plot = (x_hi-x_lo)/(N_plot-1);     
x_plot = x_lo:dx_plot:x_hi; y_plot = x_plot;
[X_plot,Y_plot] = meshgrid(x_plot,y_plot);

figure(1); clf;
Ny_FE_OTS = 20/dy_FE_OTS + 1;
idx = find( abs(X_FE_OTS - A_x*t_final) + abs(Y_FE_OTS - A_y*t_final) > 4 );
u_FE_OTS(idx) = inf;
u_plot = interp2(reshape(X_FE_OTS,Ny_FE_OTS,N+1), ...
                 reshape(Y_FE_OTS,Ny_FE_OTS,N+1), ...
                 reshape(u_FE_OTS,Ny_FE_OTS,N+1), ...
                 X_plot, Y_plot, '*nearest');
if (use_color_figures)
  mesh(X_plot, Y_plot, u_plot);
  colormap('default');
else
  % black and white
  mesh(X_plot, Y_plot, u_plot, ones(size(X_plot)));
  colormap([1 1 1; 0 0 0]);
end
axis([-10 10 -10 10 -0.1 1.1])
xlabel('x'); ylabel('y'); 
%title('Forward Euler OTS Solution');
filename = sprintf('adv_eqn_2d_FE_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(2); clf;
idx = find( sqrt( (X_FE - A_x*t_final).^2 + (Y_FE - A_y*t_final).^2 ) > 5 );
u_FE(idx) = inf;
u_plot = interp2(reshape(X_FE,N+1,N+1), ...
                 reshape(Y_FE,N+1,N+1), ...
                 reshape(u_FE,N+1,N+1), ...
                 X_plot, Y_plot, '*nearest');
if (use_color_figures)
  mesh(X_plot, Y_plot, u_plot);
  colormap('default');
else
  % black and white
  mesh(X_plot, Y_plot, u_plot, ones(size(X_plot)));
  colormap([1 1 1; 0 0 0]);
end
axis([-10 10 -10 10 -0.1 1.1])
xlabel('x'); ylabel('y'); 
%title('Forward Euler Solution');
filename = sprintf('adv_eqn_2d_FE_soln.%s', print_suffix);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display time for generating plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_total = cputime - t_start;
disp('----------------------------------------');
disp_str = sprintf('Plot Generation Time: %f', t_total);
disp(disp_str);
disp('----------------------------------------');

