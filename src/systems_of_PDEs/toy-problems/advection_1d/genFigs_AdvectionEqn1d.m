%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for a system of
% decoupled advection equations. 
%  
% Kevin T. Chu
% 2008 August
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
data_dir = 'data-advection_eqn_1d';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% set simulation parameters
debug_on = 0;
timing_on = 1;

% time integration parameters
t_init  = 0.0;
t_final = 0.5;

% grid sizes to collect data on
grid_sizes = [200 400 800 1600 3200 6400];

% allocate memory for errors
err_u_FE_OTS = zeros(size(grid_sizes));
err_v_FE_OTS = zeros(size(grid_sizes));
err_u_FE     = zeros(size(grid_sizes));
err_v_FE     = zeros(size(grid_sizes));

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
    dx = 20/N;
    dt_FE = dx/4;

    % solve advection equation using forward Euler with OTS
    disp('Forward Euler OTS');
    [u_FE_OTS, v_FE_OTS, u_exact, v_exact, x, timing_data_FE_OTS] = ...
       solveAdvectionEqn1dForwardEulerOTS(dx, ...
                                          t_init, t_final, ...
                                          debug_on, timing_on);

    % solve advection equation using forward Euler
    disp('Forward Euler');
    [u_FE, v_FE, u_exact, v_exact, x, timing_data_FE] = ...
       solveAdvectionEqn1dForwardEuler(dx, dt_FE, ...
                                       t_init, t_final, ...
                                       debug_on, timing_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ... 
         'u_FE_OTS', 'u_FE', 'u_exact', 'x', ...
         'timing_data_FE_OTS', 'timing_data_FE');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  err_u = u_FE_OTS-u_exact;
  err_u_FE_OTS(i) = norm(err_u,'inf');
  err_v = v_FE_OTS-v_exact;
  err_v_FE_OTS(i) = norm(err_v,'inf');
  err_u = u_FE-u_exact;
  err_u_FE(i) = norm(err_u,'inf');
  err_v = v_FE-v_exact;
  err_v_FE(i) = norm(err_v,'inf');

end % end loop over grid sizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate solution at low grid resolution
% for illustration purposes     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_lo_res = 200;
if (use_saved_data == 1)
    
  % case:  (use_saved_data = 1) ==> load saved solutions
           
  disp('----------------------------------------');
  disp('Loading saved data for low res solution');
  disp('----------------------------------------');
  
  filename = 'data_lo_res';
  load([data_dir, '/', filename]);
  
else

  % case:  (use_saved_data ~= 1) ==> recompute solutions
  
  disp('----------------------------------------');
  disp('Generating low res solution');
  disp('----------------------------------------');
  
  % set dx and dt
  dx = 20/N_lo_res;
  dt_FE = dx/4;

  % solve advection equation using forward Euler with OTS
  disp('Forward Euler OTS');
  [u_FE_OTS_lo_res, v_FE_OTS_lo_res, u_exact_lo_res, v_exact_lo_res, ...
   x_lo_res] = ...
     solveAdvectionEqn1dForwardEulerOTS(dx, ...
                                        t_init, t_final, ...
                                        debug_on);

  % solve advection equation using forward Euler
  disp('Forward Euler');
  [u_FE_lo_res, v_FE_lo_res, u_exact_lo_res, v_exact_lo_res, x_lo_res] = ...
     solveAdvectionEqn1dForwardEuler(dx, dt_FE, ...
                                     t_init, t_final, ...
                                     debug_on);

  % save solutions and timing data to MAT-files
  filename = 'data_lo_res';
  save([data_dir, '/', filename], ...
       'u_FE_OTS_lo_res', 'u_FE_lo_res', 'x_lo_res');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_u_FE_OTS = polyfit(log(grid_sizes),log(err_u_FE_OTS),1);
order_u_FE_OTS = -P_u_FE_OTS(1);
P_v_FE_OTS = polyfit(log(grid_sizes),log(err_v_FE_OTS),1);
order_v_FE_OTS = -P_v_FE_OTS(1);
P_u_FE = polyfit(log(grid_sizes),log(err_u_FE),1);
order_u_FE = -P_u_FE(1);
P_v_FE = polyfit(log(grid_sizes),log(err_v_FE),1);
order_v_FE = -P_v_FE(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [100 10000];
loglog(grid_sizes,err_u_FE_OTS, 'bo', ...
      'MarkerSize',14, ...
      'MarkerFaceColor','b');
hold on;
order_str = sprintf('Forward Euler (OTS)\nOrder');
text(200,5e-18,order_str);

loglog(N_plot,exp(log(N_plot)*P_u_FE(1)+P_u_FE(2)),'k');
hold on;
plot(grid_sizes,err_u_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nOrder  %1.1f', order_u_FE);
text(200,1e-5,order_str);

axis([100 10000 1e-20 1]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('advection_eqn_1d_error_u_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(2); clf;
N_plot = [100 10000];
loglog(N_plot,exp(log(N_plot)*P_v_FE_OTS(1)+P_v_FE_OTS(2)),'k');
hold on;
plot(grid_sizes,err_v_FE_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS)\nOrder = %1.1f', order_v_FE_OTS);
text(200,2e-5,order_str);

loglog(N_plot,exp(log(N_plot)*P_v_FE(1)+P_v_FE(2)),'k');
hold on;
plot(grid_sizes,err_v_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nOrder  %1.1f', order_v_FE);
text(400,0.2,order_str);

axis([100 10000 1e-6 1]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('advection_eqn_1d_error_v_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(3); clf;
plot(x_lo_res,u_FE_OTS_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
axis([-10 10 0 1.1]);
xlabel('x');
%title('Forward Euler OTS Solution u')
filename = sprintf('advection_eqn_1d_u_FE_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(4); clf;
plot(x_lo_res,v_FE_OTS_lo_res,'bo')
hold on;
plot(x,v_exact,'r')
axis([-10 10 -1 1]);
xlabel('x');
%title('Forward Euler OTS Solution v')
filename = sprintf('advection_eqn_1d_v_FE_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(5); clf;
plot(x_lo_res,u_FE_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
axis([-10 10 0 1.1]);
xlabel('x');
%title('Forward Euler Solution u')
filename = sprintf('advection_eqn_1d_u_FE_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(6); clf;
plot(x_lo_res,v_FE_lo_res,'bo')
hold on;
plot(x,v_exact,'r')
axis([-10 10 -1 1]);
xlabel('x');
%title('Forward Euler Solution v')
filename = sprintf('advection_eqn_1d_v_FE_soln.%s', print_suffix);
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

