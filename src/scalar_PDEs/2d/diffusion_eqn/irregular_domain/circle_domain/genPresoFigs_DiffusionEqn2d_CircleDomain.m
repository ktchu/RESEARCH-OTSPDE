%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for the 2d diffusion 
% equation on an circular domain for presentations.
%  
% Kevin T. Chu
% 2008 June
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add path to diffusion equation solvers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..');

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
print_format = 'png';
fig_dir = 'figures';
if ~exist(fig_dir, 'dir')
  mkdir(fig_dir);
end

% set flag for loading data from saved files (instead of recomputing solution)
use_saved_data = 1;
data_dir = 'data-diffusion_eqn_2d-circle_domain';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% set simulation parameters
debug_on  = 0;
timing_on = 1;

% physical parameters
D = 0.25;  % diffusion coefficient 

% source term type
source_term_type = 3;

% time integration parameters
% NOTE: t_init = 0.0
t_final = 0.1;

% order of extrapolant used to fill ghostcells
bc_order = 4;

% grid sizes to collect data on
grid_sizes = [50 100 200 400 800];

% allocate memory for errors
err_FE_OTS = zeros(1, length(grid_sizes));
err_FE     = zeros(1, length(grid_sizes));

% allocate memory for computation time
comp_time_FE_OTS = zeros(1, length(grid_sizes));
comp_time_FE     = zeros(1, length(grid_sizes));

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
    dx = 2/N;
    dt_FE = dx^2/4/D;

    % compute level set function that defines the domain
    % NOTE:  phi is set to be a signed distance function to make
    %        it easier to identify boundary points based on the
    %        value of phi
    N_grid = N+1;
    num_gridpts = N_grid*N_grid;
    dy = dx;
    x = -1:dx:1; y = x;
    [X,Y] = meshgrid(x,y);   % grid created with y as fastest direction
    X = reshape(X,num_gridpts,1);
    Y = reshape(Y,num_gridpts,1);

    disp('Generating level set function for boundary ...');
    theta = atan2(Y,X);
    phi = sqrt(X.^2+Y.^2) - 0.8;
    phi = reshape(phi,N_grid,N_grid);
    phi = computeDistanceFunction2d(phi,[dx,dy]);
    phi = reshape(phi,num_gridpts,1);

    % display N and bc_order
    disp('===============================');
    disp_str = sprintf('  N = %d, bc_order = %d', ...
                       N, bc_order);
    disp(disp_str);

    % solve diffusion equation using forward Euler with OTS-NIDC
    disp('---------------------------');
    disp('  Forward Euler OTS-NIDC');
    disp('---------------------------');
    [u_FE_OTS, u_exact, X, Y, timing_data_FE_OTS] = ...
       solveDiffusionEqnForwardEulerOTS2d(D, ...
                                          source_term_type, ...
                                          phi, ...
                                          dx, ...
                                          t_final, ...
                                          bc_order, ...
                                          debug_on, timing_on);

    % solve diffusion equation using forward Euler
    disp('---------------------------');
    disp('  Forward Euler');
    disp('---------------------------');
    [u_FE, u_exact, X, Y, timing_data_FE] = ...
       solveDiffusionEqnForwardEuler2d(D, ...
                                       source_term_type, ...
                                       phi, ...
                                       dx, dt_FE, ...
                                       t_final, ...
                                       bc_order, ...
                                       debug_on, timing_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ...
         'u_FE_OTS', 'u_FE', 'u_exact', ...
         'timing_data_FE_OTS', 'timing_data_FE', ...
         'phi', ...
         'X', 'Y', 'dx');
         
  end % end case:  (use_saved_data ~= 1) ==> recompute solutions
 
  % compute error
  err = u_FE_OTS-u_exact;
  err_FE_OTS(i) = norm(err,inf);
  err = u_FE-u_exact;
  err_FE(i) = norm(err,inf);

  % collect timing data
  comp_time_FE_OTS(i) = timing_data_FE_OTS(3);
  comp_time_FE(i)     = timing_data_FE(3);

end % end loop over grid sizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate solution to coarser mesh 
% for illustration purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load results generated with N = 400
N = 400;
filename = sprintf('data_%d', N);
load([data_dir, '/', filename]);

% subsample solution and convert to matrix form
N_plot = 50;
dx = 2/N_plot;
x_plot = -1:dx:1; y_plot = x_plot;
[X_plot, Y_plot] = meshgrid(x_plot, y_plot);
X_plot   = reshape(X_plot, N_plot+1, N_plot+1);
Y_plot   = reshape(Y_plot, N_plot+1, N_plot+1);
u_FE_OTS_plot = interp2(reshape(X,N+1,N+1), ...
                        reshape(Y,N+1,N+1), ...
                        reshape(u_FE_OTS,N+1,N+1), ...
                        X_plot,Y_plot,'*cubic');
err = u_FE_OTS-u_exact;
err_plot      = interp2(reshape(X,N+1,N+1), ...
                        reshape(Y,N+1,N+1), ...
                        reshape(err,N+1,N+1), ...
                        X_plot,Y_plot,'*cubic');
reshape(u_FE_OTS_plot, N_plot+1, N_plot+1);
reshape(err_plot, N_plot+1, N_plot+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_FE_OTS = polyfit(log(grid_sizes),log(err_FE_OTS),1);
P_FE = polyfit(log(grid_sizes),log(err_FE),1);
order_FE_OTS = -P_FE_OTS(1);
order_FE = -P_FE(1);


P_comp_time_FE_OTS = polyfit(log(err_FE_OTS),log(comp_time_FE_OTS),1);
comp_time_exp_FE_OTS = P_comp_time_FE_OTS(1);
P_comp_time_FE = polyfit(log(err_FE),log(comp_time_FE),1);
comp_time_exp_FE = P_comp_time_FE(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [10 1000];
loglog(N_plot, exp(log(N_plot)*P_FE_OTS(1)+P_FE_OTS(2)),'k');
hold on;
plot(grid_sizes, err_FE_OTS, 'bo', ...
     'MarkerSize', 14, ...
     'MarkerFaceColor', 'b');
order_str = sprintf('Forward Euler (OTS-NIDC)\nOrder = %1.1f', order_FE_OTS);
text(12,3e-5,order_str);

loglog(N_plot, exp(log(N_plot)*P_FE(1)+P_FE(2)), 'k');
hold on;
plot(grid_sizes, err_FE, 'rs', ...
     'MarkerSize', 14, ...
     'MarkerFaceColor', 'r');
order_str = sprintf('Forward Euler\nOrder = %1.1f', order_FE);
text(250,1e-2,order_str);

axis([10 1000 1e-8 1e-1]);
xlabel('N');
ylabel('L^\infty Error');
set(gca, 'YTick', 10.^[-8:-1]);
set(gca, 'YMinorTick', 'off');
filename = sprintf('diffusion_eqn_2d_circle_domain_error_vs_N.%s', ...
                   print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(2); clf;
err_plot = [1e-10 1e0];
loglog(err_plot, ...
  exp(log(err_plot)*P_comp_time_FE_OTS(1)+P_comp_time_FE_OTS(2)), 'k');
hold on;
loglog(err_FE_OTS, comp_time_FE_OTS, 'bo', ...
       'MarkerSize', 14, ...
       'MarkerFaceColor', 'b');
order_str = sprintf('Forward Euler (OTS-NIDC)\nOrder = %1.1f', comp_time_exp_FE_OTS);
text(3e-8,0.5,order_str);

loglog(err_plot, exp(log(err_plot)*P_comp_time_FE(1)+P_comp_time_FE(2)), 'k');
hold on;
loglog(err_FE, comp_time_FE, 'rs', ...
       'MarkerSize', 14, ...
       'MarkerFaceColor', 'r');
order_str = sprintf('Forward Euler\nOrder = %1.1f', comp_time_exp_FE);
text(3e-3,3e2,order_str);

axis([1e-8 1e0 1e-2 1e5]);
xlabel('L^\infty Error');
ylabel('Compute Time (s)');
set(gca, 'YTick', 10.^[-2:5]);
set(gca, 'YMinorTick', 'off');
filename = sprintf('diffusion_eqn_2d_circle_domain_comp_time.%s', ...
                   print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(3); clf;
surf(X_plot, Y_plot, u_FE_OTS_plot);
colormap('default');

xlabel('x'); ylabel('y'); 
axis([-1 1 -1 1 -2 2]);
view(127.5, 55);
filename = sprintf('diffusion_eqn_2d_circle_domain_soln.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(4), clf;
contour(reshape(X,N+1,N+1), ...
        reshape(Y,N+1,N+1), ...
        reshape(phi,N+1,N+1),[0 0], ...
        'LineColor', 'k', 'LineWidth', 2);
hold on
idx = find(abs(err) > 0.25*norm(err,inf));
plot(X(idx), Y(idx), 'ko', 'MarkerFaceColor', 'k');

xlabel('x'); ylabel('y'); 
filename = sprintf('diffusion_eqn_2d_circle_domain_error.%s', print_format);
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

