%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for the 2d 
% reaction-diffusion equation.
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
use_color_figures = 1;
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
use_saved_data = 1;
data_dir = 'data-RxnDiffusionEqn2d-Src';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% set simulation parameters
debug_on  = 0;
timing_on = 1;

% physical parameters
g0   = 1.2;
R_in = 1.0;

% time integration parameters
% NOTE: t_init = 0.0
t_final = 1.0;

% grid parameters
x_lo = -10;
x_hi = 10;

% grid sizes to collect data on
grid_sizes = [25 50 100 200 400];

% allocate memory for errors
err_FE_OTS = zeros(size(grid_sizes));
err_FE     = zeros(size(grid_sizes));

% allocate memory for computation time
comp_time_FE_OTS = zeros(size(grid_sizes));
comp_time_FE     = zeros(size(grid_sizes));

% start clock for timing plot generation time
t_start = cputime;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(grid_sizes)

  % get grid size
  N = grid_sizes(i);

  if (use_saved_data == 1)

    % case:  (use_saved_data = 1) ==> load saved solutions
    disp('----------------------------------------');
    disp_str = sprintf('Loading saved data for N = %d', N);
    disp(disp_str);
    disp('----------------------------------------');

    filename = sprintf('data_%d', N);
    load([data_dir, '/', filename]);

  else

    % case:  (use_saved_data ~= 1) ==> recompute solutions

    % set dx and dt
    dx = (x_hi-x_lo)/N;
    dt_FE = dx^2/4;

    disp('===============================');
    disp_str = sprintf('  N = %d', N);
    disp(disp_str);

    % solve reaction-diffusion equation using forward Euler with OTS
    disp('---------------------------');
    disp('  Forward Euler OTS')
    disp('---------------------------');
    [u_FE_OTS, u_exact_OTS, X, Y, timing_data_FE_OTS] = ...
      solveRxnDiffEqnForwardEulerOTS2d(g0, R_in, ...
                                       dx, ...
                                       t_final, ...
                                       debug_on, timing_on);

    % solve reaction-diffusion equation using forward Euler
    disp('---------------------------');
    disp('  Forward Euler')
    disp('---------------------------');
    [u_FE, u_exact, X, Y, timing_data_FE] = ...
      solveRxnDiffEqnForwardEuler2d(g0, R_in, ...
                                    dx, dt_FE, ...
                                    t_final, ...
                                    debug_on, timing_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ...
         'u_FE_OTS', 'u_FE', 'u_exact', ...
         'timing_data_FE_OTS', 'timing_data_FE', ...
         'X', 'Y');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  err = u_FE_OTS-u_exact;
  err_FE_OTS(i) = norm(err,'inf');
  err = u_FE-u_exact;
  err_FE(i) = norm(err,'inf');

  % collect timing data
  comp_time_FE_OTS(i) = timing_data_FE_OTS(2);
  comp_time_FE(i) = timing_data_FE(2);

end % end loop over grid sizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate solution to coarser mesh
% for illustration purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subsample solution and convert to matrix form
N_plot = 50;
dx = (x_hi-x_lo)/N_plot;
x_plot = x_lo:dx:x_hi; y_plot = x_plot;
[Y_plot, X_plot] = meshgrid(x_plot, y_plot);
X_plot   = reshape(X_plot, N_plot+1, N_plot+1);
Y_plot   = reshape(Y_plot, N_plot+1, N_plot+1);

u_FE_OTS_plot = interp2(reshape(X,N+1,N+1), ...
                        reshape(Y,N+1,N+1), ...
                        reshape(u_FE_OTS,N+1,N+1), ...
                        X_plot,Y_plot,'*cubic');
u_FE_plot = interp2(reshape(X,N+1,N+1), ...
                    reshape(Y,N+1,N+1), ...
                    reshape(u_FE,N+1,N+1), ...
                    X_plot,Y_plot,'*cubic');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_FE_OTS = polyfit(log(grid_sizes),log(err_FE_OTS),1);
order_FE_OTS = -P_FE_OTS(1);
P_FE = polyfit(log(grid_sizes),log(err_FE),1);
order_FE = -P_FE(1);

P_comp_time_FE_OTS = polyfit(log(err_FE_OTS(2:end)), ...
                             log(comp_time_FE_OTS(2:end)),1);
comp_time_exp_FE_OTS = P_comp_time_FE_OTS(1);
P_comp_time_FE = polyfit(log(err_FE(2:end)),log(comp_time_FE(2:end)),1);
comp_time_exp_FE = P_comp_time_FE(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_FE_OTS(1)+P_FE_OTS(2)),'k');
hold on;
plot(grid_sizes,err_FE_OTS, 'go', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','g');
order_str = sprintf('Forward Euler (OTS)\nOrder = %1.1f', order_FE_OTS);
text(18,4e-8,order_str);

loglog(N_plot,exp(log(N_plot)*P_FE(1)+P_FE(2)),'k');
hold on;
plot(grid_sizes,err_FE, 'bs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler\nOrder = %1.1f', order_FE);
text(100,3e-2,order_str);

axis([10 1000 1e-10 1e0]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('diff_eqn_2d_src_error_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(2); clf;
err_plot = [1e-10 1e0];
loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE_OTS(1)+P_comp_time_FE_OTS(2)), ...
       'k');
hold on;
loglog(err_FE_OTS, comp_time_FE_OTS, 'go', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','g');
order_str = sprintf('Forward Euler (OTS)\nSlope = %1.1f', comp_time_exp_FE_OTS);
text(5e-8,5e-2,order_str);

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE(1)+P_comp_time_FE(2)), ...
       'k');
hold on;
loglog(err_FE, comp_time_FE, 'bs', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','b');
order_str = sprintf('Forward Euler\nSlope = %1.1f', comp_time_exp_FE);
text(4e-4,3e2,order_str);

axis([1e-8 1e0 1e-4 1e4]);
xlabel('L^\infty Error');
ylabel('Compute Time');
filename = sprintf('diff_eqn_2d_src_comp_time.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(3); clf;
surf(X_plot,Y_plot,u_FE_OTS_plot);
xlabel('x'); ylabel('y'); 
%title('Forward Euler OTS Solution');
filename = sprintf('diff_eqn_2d_src_FE_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(4); clf;
surf(X_plot,Y_plot,u_FE_plot);
xlabel('x'); ylabel('y'); 
%title('Forward Euler Solution');
filename = sprintf('diff_eqn_2d_src_FE_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display time for generating plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_end = cputime;
disp('----------------------------------------');
disp_str = sprintf('Plot Generation Time: %f', t_end - t_start);
disp(disp_str);
disp('----------------------------------------');

