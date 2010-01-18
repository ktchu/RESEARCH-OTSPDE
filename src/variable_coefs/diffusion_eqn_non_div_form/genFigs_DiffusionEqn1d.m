%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for a
% 1d variable-coefficient diffusion equation with source term. 
%  
% Kevin T. Chu
% 2009 February
% Serendipity Research
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
use_saved_data = 0;
data_dir = 'data-var_coef_diffusion_eqn_1d';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% debug switch
debug_on = 0;

% grid parameters
x_lo = 0;
x_hi = 1;

% time integration parameters
% NOTE: t_init = 0.0
t_final = 0.1;

% grid sizes to collect data on
grid_sizes = [25 50 100 200 400 800]; 

% allocate memory for errors
err_FE                 = zeros(size(grid_sizes));

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

    disp('----------------------------------------');
    disp_str = sprintf('Computing solutions for N = %d', N);
    disp(disp_str);
    disp('----------------------------------------');

    % set dx and dt
    dx = (x_hi-x_lo)/N;
    dt_FE = 0.5*dx^2/9;  % 9 is max(D)

    % solve diffusion equation using FE 
    disp('FE');
    [u_FE, u_exact, x] = solveDiffusionEqn1dFE(N, dt_FE, t_final, debug_on);

    % solve diffusion equation using FE
    disp('FE_OTS_OptGrid');
    [u_FE_OTS_OptGrid, u_exact_opt_grid, x_opt_grid] = ...
      solveDiffusionEqn1dFE_OTS_OptGrid(N, t_final, debug_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ...
         'u_FE', 'u_exact', 'x', ...
         'u_FE_OTS_OptGrid', 'u_exact_opt_grid', 'x_opt_grid');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  err = u_FE-u_exact;
  err_FE(i) = norm(err,'inf');
  err = u_FE_OTS_OptGrid-u_exact_opt_grid;
  err_FE_OTS_OptGrid(i) = norm(err,'inf');

end % end loop over grid sizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate solution at low grid resolution
% for illustration purposes    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (use_saved_data == 1) 
    
  % case:  (use_saved_data = 1) ==> load saved solutions
         
  disp('----------------------------------------');
  disp('Loading saved data for low res solution');
  disp('----------------------------------------');
  
  filename = 'data_plotting';
  load([data_dir, '/', filename]);
  
else

  % case:  (use_saved_data ~= 1) ==> recompute solutions

  disp('----------------------------------------');
  disp('Generating low res solution');
  disp('----------------------------------------');

  % set dx and dt
  N_lo_res = 25;
  dx = (x_hi-x_lo)/N_lo_res;
  dt_FE = 0.5*dx^2/9;  % 9 is max(D)

  % solve diffusion equation using FE 
  disp('FE');
  [u_FE_lo_res, u_exact_lo_res, x_lo_res] = ...
    solveDiffusionEqn1dFE(N_lo_res, dt_FE, t_final, debug_on);

  % solve diffusion equation using FE-OTS-OptGrid
  disp('FE-OTS-OptGrid');
  [u_FE_OTS_OptGrid_lo_res, u_exact_lo_res_opt_grid, x_lo_res_opt_grid] = ...
    solveDiffusionEqn1dFE_OTS_OptGrid(N_lo_res, t_final, debug_on);

  disp('----------------------------------------');
  disp('Generating high res exact solution');
  disp('----------------------------------------');

  % set dx and dt
  N_hi_res = 200;
  dx = (x_hi-x_lo)/N_hi_res;
  dt_FE = 0.5*dx^2/9;  % 9 is max(D)

  % solve diffusion equation using FE
  disp('FE');
  [u_FE_hi_res, u_exact_hi_res, x_hi_res] = ...
    solveDiffusionEqn1dFE(N_hi_res, dt_FE, t_final, debug_on);

  % solve diffusion equation using FE-OTS-OptGrid
  disp('FE-OTS-OptGrid');
  [u_FE_OTS_OptGrid_hi_res, u_exact_hi_res_opt_grid, x_hi_res_opt_grid] = ...
    solveDiffusionEqn1dFE_OTS_OptGrid(N_hi_res, t_final, debug_on);

  % save solutions and timing data to MAT-files
  filename = 'data_plotting';
  save([data_dir, '/', filename], ...
       'u_FE_lo_res', 'u_FE_hi_res', ...
       'x_lo_res', 'u_exact_hi_res', 'x_hi_res', ...
       'u_FE_OTS_OptGrid_lo_res', 'u_FE_OTS_OptGrid_hi_res', ...
       'x_lo_res_opt_grid', 'u_exact_hi_res_opt_grid', 'x_hi_res_opt_grid');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_FE = polyfit(log(grid_sizes),log(err_FE),1);
order_FE = -P_FE(1);
P_FE_OTS_OptGrid = polyfit(log(grid_sizes),log(err_FE_OTS_OptGrid),1);
order_FE_OTS_OptGrid = -P_FE_OTS_OptGrid(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [10 1000];
loglog(N_plot, ...
  exp(log(N_plot)*P_FE_OTS_OptGrid(1)+P_FE_OTS_OptGrid(2)),'k');
hold on;
plot(grid_sizes,err_FE_OTS_OptGrid, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('FE-OTS-OptGrid \n(Order = %1.1f)', order_FE_OTS_OptGrid);
text(200,2e-6,order_str);

N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_FE(1)+P_FE(2)),'k');
hold on;
loglog(grid_sizes,err_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('FE (Order = %1.1f)', order_FE);
text(100,8e-2,order_str);

axis([10 1000 1e-8 1]);
set(gca, 'ytick', 10.^[-8:2:0]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('var_coef_diffusion_eqn_1d_error_vs_N.%s', print_suffix);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(2); clf;
plot(x_lo_res,u_FE_lo_res,'bo')
hold on;
plot(x_hi_res,u_exact_hi_res,'r')
axis([x_lo x_hi 0 4]);
xlabel('x');
%title('FE Solution')
filename = sprintf('var_coef_diffusion_eqn_1d_FE_soln.%s', print_suffix);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(3); clf;
plot(x_lo_res_opt_grid,u_FE_OTS_OptGrid_lo_res,'bo')
hold on;
plot(x_hi_res,u_exact_hi_res,'r')
axis([x_lo x_hi 0 4]);
xlabel('x');
%title('FE OTS-OptGrid Solution')
filename = sprintf('var_coef_diffusion_eqn_1d_FE_OTS_OptGridsoln.%s', ...
  print_suffix);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display time for generating plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_end = cputime;
disp('----------------------------------------');
disp_str = sprintf('Plot Generation Time: %f', t_end-t_start);
disp(disp_str);
disp('----------------------------------------');

