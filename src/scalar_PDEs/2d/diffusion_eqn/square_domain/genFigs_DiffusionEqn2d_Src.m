%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for the 2d diffusion 
% equation without a src term.
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
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontName','Helvetica')
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultTextFontSize',16,'DefaultTextFontName','Helvetica')

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
data_dir = 'data-diffusion_eqn_2d_src';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% set simulation parameters
debug_on  = 0;
timing_on = 1;

% src term type
src_term_type = 3;

% physical parameters
D = 0.25;  % diffusion coefficient 

% time integration parameters
% NOTE: t_init = 0.0
t_final = 0.1;

% grid sizes to collect data on
grid_sizes = [25 50 100 200 400 800];

% allocate memory for errors
err_FE_OTS = zeros(size(grid_sizes));
err_FE     = zeros(size(grid_sizes));
err_CN     = zeros(size(grid_sizes));

% allocate memory for computation time
comp_time_FE_OTS = zeros(size(grid_sizes));
comp_time_FE     = zeros(size(grid_sizes));
comp_time_CN     = zeros(size(grid_sizes));

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
    dx = 1/N;
    dt_FE = dx^2/4/D;
    dt_CN = dx/16/D;

    disp('===============================');
    disp_str = sprintf('  N = %d', N);
    disp(disp_str);

    % solve diffusion equation using forward Euler with OTS
    if (N < 800)   % skip N = 800 case -- takes ~ 25-30h to compute
      disp('---------------------------');
      disp('  Forward Euler OTS')
      disp('---------------------------');
      [u_FE_OTS, u_exact, X, Y, timing_data_FE_OTS] = ...
         solveDiffusionEqnForwardEulerOTS2d(D, ...
                                            src_term_type, ...
                                            dx, ...
                                            t_final, ...
                                            debug_on, timing_on);
    end

    % solve diffusion equation using forward Euler
    disp('---------------------------');
    disp('  Forward Euler')
    disp('---------------------------');
    [u_FE, u_exact, X, Y, timing_data_FE] = ...
       solveDiffusionEqnForwardEuler2d(D, ...
                                       src_term_type, ...
                                       dx, dt_FE, ...
                                       t_final, ...
                                       debug_on, timing_on);

    % solve diffusion equation using Crank-Nicholson
    disp('---------------------------');
    disp('  Crank-Nicholson')
    disp('---------------------------');
    [u_CN, u_exact, X, Y, timing_data_CN] = ...
       solveDiffusionEqnCrankNicholson2d(D, ...
                                         src_term_type, ...
                                         dx, dt_CN, ...
                                         t_final, ...
                                         debug_on, timing_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ...
         'u_FE_OTS', 'u_FE', 'u_CN', 'u_exact', ...
         'timing_data_FE_OTS', 'timing_data_FE', 'timing_data_CN', ...
         'X', 'Y');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  if (N < 800)   % skip N = 800 case -- takes ~ 25-30h to compute
    err = u_FE_OTS-u_exact;
    err_FE_OTS(i) = norm(err,'inf');
  end
  err = u_FE-u_exact;
  err_FE(i) = norm(err,'inf');
  err = u_CN-u_exact;
  err_CN(i) = norm(err,'inf');

  % collect timing data
  if (N < 800)   % skip N = 800 case -- takes ~ 25-30h to compute
    comp_time_FE_OTS(i) = timing_data_FE_OTS(2);
  end
  comp_time_FE(i) = timing_data_FE(2);
  comp_time_CN(i) = timing_data_CN(2);

end % end loop over grid sizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate solution to coarser mesh
% for illustration purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subsample solution and convert to matrix form
N_plot = 50;
dx = 1/N_plot;
x_plot = 0:dx:1; y_plot = x_plot;
[X_plot, Y_plot] = meshgrid(x_plot, y_plot);

% FE_OTS requires special handling because it doesn't have solution 
% with N = 800
N_FE_OTS = 401;
dx_FE_OTS = 1/(N_FE_OTS-1);
x_FE_OTS = 0:dx_FE_OTS:1; y_FE_OTS = x_FE_OTS;
[X_FE_OTS, Y_FE_OTS] = meshgrid(x_FE_OTS, y_FE_OTS);
u_FE_OTS_plot = interp2(reshape(X_FE_OTS,N_FE_OTS,N_FE_OTS), ...
                        reshape(Y_FE_OTS,N_FE_OTS,N_FE_OTS), ...
                        reshape(u_FE_OTS,N_FE_OTS,N_FE_OTS), ...
                        X_plot,Y_plot,'*cubic');
u_FE_plot = interp2(reshape(X,N+1,N+1), ...
                    reshape(Y,N+1,N+1), ...
                    reshape(u_FE,N+1,N+1), ...
                    X_plot,Y_plot,'*cubic');
u_CN_plot = interp2(reshape(X,N+1,N+1), ...
                    reshape(Y,N+1,N+1), ...
                    reshape(u_CN,N+1,N+1), ...
                    X_plot,Y_plot,'*cubic');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_FE_OTS = polyfit(log(grid_sizes(1:end-1)),log(err_FE_OTS(1:end-1)),1);
order_FE_OTS = -P_FE_OTS(1);
P_FE = polyfit(log(grid_sizes),log(err_FE),1);
order_FE = -P_FE(1);
P_CN = polyfit(log(grid_sizes),log(err_CN),1);
order_CN = -P_CN(1);

P_comp_time_FE_OTS = polyfit(log(err_FE_OTS(1:end-1)), ...
                             log(comp_time_FE_OTS(1:end-1)),1);
comp_time_exp_FE_OTS = P_comp_time_FE_OTS(1);
P_comp_time_FE = polyfit(log(err_FE),log(comp_time_FE),1);
comp_time_exp_FE = P_comp_time_FE(1);
P_comp_time_CN = polyfit(log(err_CN),log(comp_time_CN),1);
comp_time_exp_CN = P_comp_time_CN(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_FE_OTS(1)+P_FE_OTS(2)),'k');
hold on;
plot(grid_sizes,err_FE_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS-NIDC)\nOrder = %1.1f', order_FE_OTS);
text(20,1e-8,order_str);

loglog(N_plot,exp(log(N_plot)*P_FE(1)+P_FE(2)),'k');
hold on;
plot(grid_sizes,err_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nOrder = %1.1f', order_FE);
text(40,1e-1,order_str);

loglog(N_plot,exp(log(N_plot)*P_CN(1)+P_CN(2)),'k');
hold on;
plot(grid_sizes,err_CN, 'md', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','m');
order_str = sprintf('Crank-Nicholson\nOrder = %1.1f', order_CN);
text(150,1e-2,order_str);

axis([10 1000 1e-10 1e0]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('diffusion_eqn_2d_src_error_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


f = figure(2); clf;
err_plot = [1e-10 1e0];
loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE_OTS(1)+P_comp_time_FE_OTS(2)), ...
       'k');
hold on;
loglog(err_FE_OTS, comp_time_FE_OTS, 'bo', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS-NIDC)\nSlope = %1.1f', comp_time_exp_FE_OTS);
text(3e-7,3e-3,order_str);

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE(1)+P_comp_time_FE(2)), ...
       'k');
hold on;
loglog(err_FE, comp_time_FE, 'rs', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nSlope = %1.1f', comp_time_exp_FE);
text(5e-4,5e2,order_str);

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_CN(1)+P_comp_time_CN(2)), ...
       'k');
hold on;
loglog(err_CN, comp_time_CN, 'md', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','m');
order_str = sprintf('Crank-Nicholson\nSlope = %1.1f', comp_time_exp_CN);
text(4e-8,1.5e3,order_str);
annotation(f, 'arrow', [0.36 0.45], [0.82 0.82], ...
           'linewidth', 3, 'headstyle', 'plain');

axis([1e-8 1e0 1e-4 1e4]);
xlabel('L^\infty Error');
ylabel('Compute Time (s)');
filename = sprintf('diffusion_eqn_2d_src_comp_time.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(3); clf;
surf(x_plot,y_plot,u_FE_OTS_plot);
xlabel('x'); ylabel('y'); 
%title('Forward Euler OTS Solution');
filename = sprintf('diffusion_eqn_2d_src_FE_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(4); clf;
surf(x_plot,y_plot,u_FE_plot);
xlabel('x'); ylabel('y'); 
%title('Forward Euler Solution');
filename = sprintf('diffusion_eqn_2d_src_FE_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(5); clf;
surf(x_plot,y_plot,u_CN_plot);
xlabel('x'); ylabel('y'); 
%title('Crank-Nicholson Solution');
filename = sprintf('diffusion_eqn_2d_src_CN_soln.%s', print_suffix);
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

