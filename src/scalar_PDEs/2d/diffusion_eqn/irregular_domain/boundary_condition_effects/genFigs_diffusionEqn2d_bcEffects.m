%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures to compare the accuracy achieved 
% using various orders of boundary conditions for the 2d diffusion equation 
% on an irregular domain.
%  
% Kevin T. Chu
% 2008 February
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
use_color_mesh_figures = 1;

% domain geometry (valid choices: 'circle', 'starfish')
domain_geometry = 'circle';

% set flag for loading data from saved files (instead of recomputing solution)
use_saved_data = 1;
if strcmp(domain_geometry, 'circle')
  data_dir = 'data-diffusion_eqn_2d-circle_domain';
elseif strcmp(domain_geometry, 'starfish')
  data_dir = 'data-diffusion_eqn_2d-starfish_domain';
else
  error('Invalid domain geometry')
end
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

% grid sizes to collect data on
grid_sizes = [50 100 200 400];
grid_sizes = [50 100];

% allocate memory for errors
err_FE_OTS = zeros(length(grid_sizes), 3);
err_FE     = zeros(length(grid_sizes), 1);

% allocate memory for computation time
comp_time_FE_OTS = zeros(length(grid_sizes), 3);
comp_time_FE     = zeros(length(grid_sizes), 1);

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
    dt_FE = dx^2/8/D;

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
    if strcmp(domain_geometry, 'circle')
      phi = sqrt(X.^2+Y.^2) - 0.8;                        % circle
    elseif strcmp(domain_geometry, 'starfish')
      phi = sqrt(X.^2+Y.^2) - 0.6*(1+0.5*sin(5*theta));   % starfish
    else
      error('Invalid domain geometry')
    end
    phi = reshape(phi,N_grid,N_grid);
    phi = computeDistanceFunction2d(phi,[dx,dy]);
    phi = reshape(phi,num_gridpts,1);

    % compute solution with bc_order = 2
    bc_order = 2;
    disp('===============================');
    disp_str = sprintf('  N = %d, bc_order = %d', ...
                       N, bc_order);
    disp(disp_str);

    % solve diffusion equation using forward Euler with OTS
    disp('---------------------------');
    disp_str = sprintf('  Forward Euler OTS, BC Order = %d', bc_order);
    disp(disp_str);
    disp('---------------------------');
    [u_FE_OTS_bc_order_2, u_exact, X, Y, timing_data_FE_OTS_bc_order_2] = ...
       solveDiffusionEqnForwardEulerOTS2d(D, ...
                                          source_term_type, ...
                                          phi, ...
                                          dx, ...
                                          t_final, ...
                                          bc_order, ...
                                          debug_on, timing_on);

    % compute solution with bc_order = 3
    bc_order = 3;
    disp('===============================');
    disp_str = sprintf('  N = %d, bc_order = %d', ...
                       N, bc_order);
    disp(disp_str);

    % solve diffusion equation using forward Euler with OTS
    disp('---------------------------');
    disp_str = sprintf('  Forward Euler OTS, BC Order = %d', bc_order);
    disp(disp_str);
    disp('---------------------------');
    [u_FE_OTS_bc_order_3, u_exact, X, Y, timing_data_FE_OTS_bc_order_3] = ...
       solveDiffusionEqnForwardEulerOTS2d(D, ...
                                          source_term_type, ...
                                          phi, ...
                                          dx, ...
                                          t_final, ...
                                          bc_order, ...
                                          debug_on, timing_on);

    % compute solution with bc_order = 4
    bc_order = 4;
    disp('===============================');
    disp_str = sprintf('  N = %d, bc_order = %d', ...
                       N, bc_order);
    disp(disp_str);

    % solve diffusion equation using forward Euler with OTS
    disp('---------------------------');
    disp_str = sprintf('  Forward Euler OTS, BC Order = %d', bc_order);
    disp(disp_str);
    disp('---------------------------');
    [u_FE_OTS_bc_order_4, u_exact, X, Y, timing_data_FE_OTS_bc_order_4] = ...
       solveDiffusionEqnForwardEulerOTS2d(D, ...
                                          source_term_type, ...
                                          phi, ...
                                          dx, ...
                                          t_final, ...
                                          bc_order, ...
                                          debug_on, timing_on);

    % solve diffusion equation using forward Euler
    bc_order = 4;
    disp('---------------------------');
    disp_str = sprintf('  Forward Euler, BC Order = %d', bc_order);
    disp(disp_str);
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
         'u_FE_OTS_bc_order_2', 'u_FE_OTS_bc_order_3', ...
         'u_FE_OTS_bc_order_4', 'u_FE', 'u_exact', ...
         'timing_data_FE_OTS_bc_order_2', 'timing_data_FE_OTS_bc_order_3', ...
         'timing_data_FE_OTS_bc_order_4', 'timing_data_FE', ...
         'phi', ...
         'X', 'Y', 'dx');
         
  end % end case:  (use_saved_data ~= 1) ==> recompute solutions
 
  % compute error
  err = u_FE_OTS_bc_order_2-u_exact;
  err_FE_OTS(i, 1) = norm(err,inf);
  err = u_FE_OTS_bc_order_3-u_exact;
  err_FE_OTS(i, 2) = norm(err,inf);
  err = u_FE_OTS_bc_order_4-u_exact;
  err_FE_OTS(i, 3) = norm(err,inf);
  err = u_FE-u_exact;
  err_FE(i) = norm(err,inf);

  % collect timing data
  comp_time_FE_OTS(i, 1) = timing_data_FE_OTS_bc_order_2(3);
  comp_time_FE_OTS(i, 2) = timing_data_FE_OTS_bc_order_3(3);
  comp_time_FE_OTS(i, 3) = timing_data_FE_OTS_bc_order_4(3);
  comp_time_FE(i)        = timing_data_FE(3);

end % end loop over grid sizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate solution to coarser mesh 
% for illustration purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subsample solution and convert to matrix form
N_plot = 50;
dx = 2/N_plot;
x_plot = -1:dx:1; y_plot = x_plot;
[X_plot, Y_plot] = meshgrid(x_plot, y_plot);
X_plot   = reshape(X_plot, N_plot+1, N_plot+1);
Y_plot   = reshape(Y_plot, N_plot+1, N_plot+1);
u_FE_OTS_plot = interp2(reshape(X,N+1,N+1), ...
                        reshape(Y,N+1,N+1), ...
                        reshape(u_FE_OTS_bc_order_4,N+1,N+1), ...
                        X_plot,Y_plot,'*cubic');
err = u_FE_OTS_bc_order_4-u_exact;
err_plot      = interp2(reshape(X,N+1,N+1), ...
                        reshape(Y,N+1,N+1), ...
                        reshape(err,N+1,N+1), ...
                        X_plot,Y_plot,'*cubic');
reshape(u_FE_OTS_plot, N_plot+1, N_plot+1);
reshape(err_plot, N_plot+1, N_plot+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_FE_OTS_bc_order_2 = polyfit(log(grid_sizes),log(err_FE_OTS(:,1)'),1);
P_FE_OTS_bc_order_3 = polyfit(log(grid_sizes),log(err_FE_OTS(:,2)'),1);
P_FE_OTS_bc_order_4 = polyfit(log(grid_sizes),log(err_FE_OTS(:,3)'),1);
P_FE_bc_order_4 = polyfit(log(grid_sizes),log(err_FE'),1);
order_FE_OTS_bc_order_2 = -P_FE_OTS_bc_order_2(1);
order_FE_OTS_bc_order_3 = -P_FE_OTS_bc_order_3(1);
order_FE_OTS_bc_order_4 = -P_FE_OTS_bc_order_4(1);
order_FE_bc_order_4 = -P_FE_bc_order_4(1);


P_comp_time_FE_OTS_bc_order_2 = ...
   polyfit(log(err_FE_OTS(:,1)),log(comp_time_FE_OTS(:,1)),1);
comp_time_exp_FE_OTS_bc_order_2 = P_comp_time_FE_OTS_bc_order_2(1);
P_comp_time_FE_OTS_bc_order_3 = ...
   polyfit(log(err_FE_OTS(:,2)),log(comp_time_FE_OTS(:,2)),1);
comp_time_exp_FE_OTS_bc_order_3 = P_comp_time_FE_OTS_bc_order_3(1);
P_comp_time_FE_OTS_bc_order_4 = ...
   polyfit(log(err_FE_OTS(:,3)),log(comp_time_FE_OTS(:,3)),1);
comp_time_exp_FE_OTS_bc_order_4 = P_comp_time_FE_OTS_bc_order_4(1);
P_comp_time_FE_bc_order_4 = ...
   polyfit(log(err_FE),log(comp_time_FE),1);
comp_time_exp_FE = P_comp_time_FE_bc_order_4(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [10 1000];
loglog(N_plot, ...
       exp(log(N_plot)*P_FE_OTS_bc_order_4(1)+P_FE_OTS_bc_order_4(2)),'k');
hold on;
plot(grid_sizes, err_FE_OTS(:,3), 'go', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','g');
order_str = sprintf('Forward Euler OTS\n4th-order BC\nOrder = %1.1f', ...
                     order_FE_OTS_bc_order_4);
text(12,1e-3,order_str);

loglog(N_plot, ...
       exp(log(N_plot)*P_FE_OTS_bc_order_3(1)+P_FE_OTS_bc_order_3(2)),'k');
hold on;
plot(grid_sizes, err_FE_OTS(:,2), 'rd', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler OTS\n3rd-order BC\nOrder = %1.1f', order_FE_OTS_bc_order_3);
text(30,1e-5,order_str);

loglog(N_plot, ...
       exp(log(N_plot)*P_FE_OTS_bc_order_2(1)+P_FE_OTS_bc_order_2(2)),'k');
hold on;
plot(grid_sizes, err_FE_OTS(:,1), 'bs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler OTS\n2nd-order BC\nOrder = %1.1f', ...
                    order_FE_OTS_bc_order_2);
text(70,2e-1,order_str);

loglog(N_plot, ...
       exp(log(N_plot)*P_FE_bc_order_4(1)+P_FE_bc_order_4(2)),'k');
hold on;
plot(grid_sizes, err_FE, 'mx', ...
     'MarkerSize',14);
order_str = sprintf('Forward Euler\n4th-order BC\nOrder = %1.1f', ...
                    order_FE_bc_order_4);
text(250,1e-2,order_str);

axis([10 1000 1e-8 1e0]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('diffusion_eqn_2d_irreg_domain_error_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(2); clf;
err_plot = [1e-10 1e0];
loglog(err_plot, ...
  exp(log(err_plot)*P_comp_time_FE_OTS_bc_order_4(1) ...
     +P_comp_time_FE_OTS_bc_order_4(2)), ...
  'k');
hold on;
loglog(err_FE_OTS(:,3), comp_time_FE_OTS(:,3), 'go', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','g');
order_str = sprintf('Forward Euler OTS\n4th-order BC\nOrder = %1.1f', ...
                     comp_time_exp_FE_OTS_bc_order_4);
text(5e-8,5e-2,order_str);

loglog(err_plot, ...
  exp(log(err_plot)*P_comp_time_FE_OTS_bc_order_3(1) ...
     +P_comp_time_FE_OTS_bc_order_3(2)), ...
  'k');
hold on;
loglog(err_FE_OTS(:,2), comp_time_FE_OTS(:,2), 'rd', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','r');
order_str = sprintf('Forward Euler OTS\n3rd-order BC\nOrder = %1.1f', ...
                     comp_time_exp_FE_OTS_bc_order_3);
text(5e-7,1e0,order_str);

loglog(err_plot, ...
  exp(log(err_plot)*P_comp_time_FE_OTS_bc_order_2(1) ...
     +P_comp_time_FE_OTS_bc_order_2(2)), ...
  'k');
hold on;
loglog(err_FE_OTS(:,1), comp_time_FE_OTS(:,3), 'bs', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','b');
order_str = sprintf('Forward Euler OTS\n2nd-order BC\nOrder = %1.1f', ...
                     comp_time_exp_FE_OTS_bc_order_2);
text(5e-6,100,order_str);

loglog(err_plot, ...
  exp(log(err_plot)*P_comp_time_FE_bc_order_4(1) ...
     +P_comp_time_FE_bc_order_4(2)), ...
  'k');
hold on;
loglog(err_FE, comp_time_FE, 'mx', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','m');
order_str = sprintf('Forward Euler\n4th-order BC\nOrder = %1.1f', ...
                    comp_time_exp_FE);
text(4e-4,3e2,order_str);

axis([1e-8 1e0 1e-4 1e4]);
xlabel('L^\infty Error');
ylabel('Compute Time (s)');
filename = sprintf('diffusion_eqn_2d_src_comp_time.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(3); clf;
if (use_color_mesh_figures)
  mesh(X_plot, Y_plot, u_FE_OTS_plot);
  colormap('default');
else
  % black and white
  mesh(X_plot, Y_plot, u_FE_OTS_plot, ones(size(X_plot)));
  colormap([1 1 1; 0 0 0]);
end

xlabel('x'); ylabel('y'); 
axis([-1 1 -1 1 -2 2]);
view(127.5, 55);
filename = sprintf('diffusion_eqn_2d_irreg_domain_soln.%s', print_suffix);
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

