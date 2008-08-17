%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for a
% 1d variable-coefficient wave equation with source term. 
%  
% Kevin T. Chu
% 2008 August
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
use_saved_data = 1;
data_dir = 'data-var_coef_wave_eqn_1d';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% debug switch
debug_on = 0;

% grid parameters
x_lo = -1;
x_hi = 1;

% time integration parameters
% NOTE: t_init = 0.0
t_final = 1.0;

% grid sizes to collect data on
grid_sizes = [100 200 400 800 1600 3200 6400];

% allocate memory for errors
err_KPY_Transformed_OTS = zeros(size(grid_sizes));
err_KPY_Transformed     = zeros(size(grid_sizes));
err_KPY                 = zeros(size(grid_sizes));

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
    dt_KPY = dx/4;

    % solve wave equation using KPY on a transformed domain with OTS
    disp('KPY-Transformed-OTS');
    [u_KPY_Transformed_OTS, u_exact_Transformed, x_Transformed] = ...
      solveWaveEqn1dKPY_Transformed_OTS(N, t_final, debug_on);

    % solve wave equation using KPY on a transformed domain without OTS
    disp('KPY-Transformed');
    [u_KPY_Transformed, u_exact_Transformed, x_Transformed] = ...
      solveWaveEqn1dKPY_Transformed(N, dt_KPY, t_final, debug_on);

    % solve wave equation using KPY on a untransformed domain 
    disp('KPY');
    [u_KPY, u_exact, x] = solveWaveEqn1dKPY(N, dt_KPY, t_final, debug_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ...
         'u_KPY_Transformed_OTS', 'u_KPY_Transformed', 'u_KPY', ...
         'u_exact_Transformed', 'u_exact', ...
         'x_Transformed', 'x');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  err = u_KPY_Transformed_OTS-u_exact_Transformed;
  err_KPY_Transformed_OTS(i) = norm(err,'inf');
  err = u_KPY_Transformed-u_exact_Transformed;
  err_KPY_Transformed(i) = norm(err,'inf');
  err = u_KPY-u_exact;
  err_KPY(i) = norm(err,'inf');

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
  N_lo_res = 100;
  dx = (x_hi-x_lo)/N_lo_res;
  dt_KPY = dx/4;

  % solve wave equation using KPY on transformed domain with OTS
  disp('KPY-Transformed-OTS');
  [u_KPY_Transformed_OTS_lo_res, u_exact_Transformed_lo_res, ...
   x_Transformed_lo_res] = ...
    solveWaveEqn1dKPY_Transformed_OTS(N_lo_res, t_final, debug_on);

  % solve wave equation using KPY on transformed domain without OTS
  disp('KPY-Transformed');
  [u_KPY_Transformed_lo_res, u_exact_Transformed_lo_res, ...
   x_Transformed_lo_res] = ...
    solveWaveEqn1dKPY_Transformed(N_lo_res, dt_KPY, t_final, debug_on);

  % solve wave equation using KPY on untransformed domain
  disp('KPY');
  [u_KPY_lo_res, u_exact_lo_res, x_lo_res] = ...
    solveWaveEqn1dKPY(N_lo_res, dt_KPY, t_final, debug_on);

  disp('----------------------------------------');
  disp('Generating high res exact solution');
  disp('----------------------------------------');

  % set dx and dt
  N_hi_res = 400;
  dx = (x_hi-x_lo)/N_hi_res;
  dt_KPY = dx/4;

  % solve wave equation using KPY on untransformed domain
  [u_KPY_hi_res, u_exact_hi_res, x_hi_res] = ...
    solveWaveEqn1dKPY(N_hi_res, dt_KPY, t_final, debug_on);

  % save solutions and timing data to MAT-files
  filename = 'data_plotting';
  save([data_dir, '/', filename], ...
       'u_KPY_Transformed_OTS_lo_res', 'u_KPY_Transformed_lo_res', ...
       'u_KPY_lo_res', ...
       'x_lo_res', 'x_Transformed_lo_res', ...
       'u_exact_hi_res', 'x_hi_res');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_KPY_Transformed_OTS = ...
  polyfit(log(grid_sizes(1:end-2)),log(err_KPY_Transformed_OTS(1:end-2)),1);
order_KPY_Transformed_OTS = -P_KPY_Transformed_OTS(1);
P_KPY_Transformed = ...
  polyfit(log(grid_sizes(1:end-2)),log(err_KPY_Transformed(1:end-2)),1);
order_KPY_Transformed = -P_KPY_Transformed(1);
P_KPY = polyfit(log(grid_sizes),log(err_KPY),1);
order_KPY = -P_KPY(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [100 10000];
loglog(N_plot, ...
  exp(log(N_plot)*P_KPY_Transformed_OTS(1)+P_KPY_Transformed_OTS(2)),'k');
hold on;
plot(grid_sizes,err_KPY_Transformed_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('KPY-Transformed-OTS\nOrder = %1.1f', ...
  order_KPY_Transformed_OTS);
text(200,2e-8,order_str);

N_plot = [100 10000];
loglog(N_plot, ...
  exp(log(N_plot)*P_KPY_Transformed(1)+P_KPY_Transformed(2)),'k');
hold on;
plot(grid_sizes,err_KPY_Transformed, 'g^', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','g');

N_plot = [100 10000];
loglog(N_plot,exp(log(N_plot)*P_KPY(1)+P_KPY(2)),'k');
hold on;
plot(grid_sizes,err_KPY, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('KPY, KPY-Transformed\nOrder = %1.1f', order_KPY);
text(800,5e-3,order_str);

axis([100 10000 1e-9 1e-1]);
set(gca, 'ytick', 10.^[-9:2:-1]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('var_coef_wave_eqn_1d_error_vs_N.%s', print_suffix);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(2); clf;
plot(x_Transformed_lo_res,u_KPY_Transformed_OTS_lo_res,'bo')
hold on;
plot(x_hi_res,u_exact_hi_res,'r')
axis([x_lo x_hi -5 5]);
xlabel('x');
%title('KPY-Transformed OTS Solution')
filename = sprintf('var_coef_wave_eqn_1d_KPY_Transformed_OTS_soln.%s', ...
                   print_suffix);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(3); clf;
plot(x_Transformed_lo_res,u_KPY_Transformed_lo_res,'bo')
hold on;
plot(x_hi_res,u_exact_hi_res,'r')
axis([x_lo x_hi -5 5]);
xlabel('x');
%title('KPY-Transformed Solution')
filename = sprintf('var_coef_wave_eqn_1d_KPY_Transformed_soln.%s', ...
                   print_suffix);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(4); clf;
plot(x_lo_res,u_KPY_lo_res,'bo')
hold on;
plot(x_hi_res,u_exact_hi_res,'r')
axis([x_lo x_hi -5 5]);
xlabel('x');
%title('KPY Solution')
filename = sprintf('var_coef_wave_eqn_1d_KPY_soln.%s', print_suffix);
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

