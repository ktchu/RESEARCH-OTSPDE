%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for a system of
% decoupled diffusion equation. 
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
use_saved_data = 1;
data_dir = 'data-system_of_diffusion_eqns_1d';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% set simulation parameters
debug_on = 0;
timing_on = 1;

% time integration parameters
t_init  = 0.0;
t_final = 0.1;

% grid sizes to collect data on
grid_sizes = [25 50 100 200 400 800 1600];

% allocate memory for errors
err_u_FE_OTS = zeros(size(grid_sizes));
err_v_FE_OTS = zeros(size(grid_sizes));
err_FE_OTS   = zeros(size(grid_sizes));
err_u_FE     = zeros(size(grid_sizes));
err_v_FE     = zeros(size(grid_sizes));
err_FE       = zeros(size(grid_sizes));

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
    dx = 2/N;
    dt_FE = dx^2/3;

    % solve diffusion equation using forward Euler with OTS
    disp('Forward Euler OTS');
    [u_FE_OTS, v_FE_OTS, u_exact, v_exact, x, timing_data_FE_OTS] = ...
       solveDiffusionEqn1dForwardEulerOTS(dx, ...
                                          t_init, t_final, ...
                                          debug_on, timing_on);

    % solve diffusion equation using forward Euler
    disp('Forward Euler');
    [u_FE, v_FE, u_exact, v_exact, x, timing_data_FE] = ...
       solveDiffusionEqn1dForwardEuler(dx, dt_FE, ...
                                       t_init, t_final, ...
                                       debug_on, timing_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ... 
         'u_FE_OTS', 'u_FE', 'u_exact', ...
         'v_FE_OTS', 'v_FE', 'v_exact', ...
         'x', 'timing_data_FE_OTS', 'timing_data_FE');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  err_u = u_FE_OTS-u_exact;
  err_u_FE_OTS(i) = norm(err_u,'inf');
  err_v = v_FE_OTS-v_exact;
  err_v_FE_OTS(i) = norm(err_v,'inf');
  err_FE_OTS(i) = max(err_u_FE_OTS(i), err_v_FE_OTS(i));
  err_u = u_FE-u_exact;
  err_u_FE(i) = norm(err_u,'inf');
  err_v = v_FE-v_exact;
  err_v_FE(i) = norm(err_v,'inf');
  err_FE(i) = max(err_u_FE(i), err_v_FE(i));

  % collect timing data
  comp_time_FE_OTS(i) = timing_data_FE_OTS;
  comp_time_FE(i) = timing_data_FE;

end % end loop over grid sizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate solution at low grid resolution
% for illustration purposes     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_lo_res = 50;
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
  dx = 2/N_lo_res;
  dt_FE = dx^2/3;

  % solve diffusion equation using forward Euler with OTS
  disp('Forward Euler OTS');
  [u_FE_OTS_lo_res, v_FE_OTS_lo_res, ...
   u_exact_lo_res, v_exact_lo_res, ...
   x_lo_res] = ...
     solveDiffusionEqn1dForwardEulerOTS(dx, ...
                                        t_init, t_final, ...
                                        debug_on);

  % solve diffusion equation using forward Euler
  disp('Forward Euler');
  [u_FE_lo_res, v_FE_lo_res, u_exact_lo_res, v_exact_lo_res, x_lo_res] = ...
     solveDiffusionEqn1dForwardEuler(dx, dt_FE, ...
                                     t_init, t_final, ...
                                     debug_on);

  % save solutions and timing data to MAT-files
  filename = 'data_lo_res';
  save([data_dir, '/', filename], ...
       'u_FE_OTS_lo_res', 'v_FE_OTS_lo_res', 'u_FE_lo_res', 'v_FE_lo_res', ...
       'x_lo_res');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_u_FE_OTS = polyfit(log(grid_sizes(2:end)),log(err_u_FE_OTS(2:end)),1);
order_u_FE_OTS = -P_u_FE_OTS(1);
P_v_FE_OTS = polyfit(log(grid_sizes(2:end)),log(err_v_FE_OTS(2:end)),1);
order_v_FE_OTS = -P_v_FE_OTS(1);
P_u_FE = polyfit(log(grid_sizes(2:end)),log(err_u_FE(2:end)),1);
order_u_FE = -P_u_FE(1);
P_v_FE = polyfit(log(grid_sizes(2:end)),log(err_v_FE(2:end)),1);
order_v_FE = -P_v_FE(1);

P_comp_time_FE_OTS = polyfit(log(err_FE_OTS(3:end)), ...
                             log(comp_time_FE_OTS(3:end)),1);
comp_time_exp_FE_OTS = P_comp_time_FE_OTS(1);
P_comp_time_FE = polyfit(log(err_FE(3:end)),log(comp_time_FE(3:end)),1);
comp_time_exp_FE = P_comp_time_FE(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_u_FE_OTS(1)+P_u_FE_OTS(2)),'k');
hold on;
plot(grid_sizes,err_u_FE_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS)\nOrder = %1.1f', order_u_FE_OTS);
text(20,1.e-8,order_str);

loglog(N_plot,exp(log(N_plot)*P_u_FE(1)+P_u_FE(2)),'k');
hold on;
plot(grid_sizes,err_u_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nOrder = %1.1f', order_u_FE);
text(100,1e-2,order_str);

axis([10 1000 1e-10 1]);
xlabel('N');
ylabel('L^\infty Error');
set(gca,'ytick',10.^[-10:2:0]);
filename = sprintf('system_of_diffusion_eqns_1d_error_u_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(2); clf;
N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_v_FE_OTS(1)+P_v_FE_OTS(2)),'k');
hold on;
plot(grid_sizes,err_v_FE_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS)\nOrder = %1.1f', order_v_FE_OTS);
text(20,1.e-8,order_str);

loglog(N_plot,exp(log(N_plot)*P_v_FE(1)+P_v_FE(2)),'k');
hold on;
plot(grid_sizes,err_v_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nOrder = %1.1f', order_v_FE);
text(100,1e-2,order_str);

axis([10 1000 1e-10 1]);
xlabel('N');
ylabel('L^\infty Error');
set(gca,'ytick',10.^[-10:2:0]);
filename = sprintf('system_of_diffusion_eqns_1d_error_v_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(3); clf;
err_plot = [1e-12 1e0];
loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE_OTS(1)+P_comp_time_FE_OTS(2)), ...
       'k');
hold on;
loglog(err_FE_OTS(3:end), comp_time_FE_OTS(3:end), 'bo', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS)\nSlope = %1.1f', comp_time_exp_FE_OTS);
text(4e-10,4e-3,order_str);

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE(1)+P_comp_time_FE(2)), ...
       'k');
hold on;
loglog(err_FE(3:end), comp_time_FE(3:end), 'rs', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nSlope = %1.1f', comp_time_exp_FE);
text(1e-7,2,order_str);

axis([1e-10 1e-2 1e-3 1e2]);
xlabel('L^\infty Error');
ylabel('Compute Time (s)');
set(gca,'ytick',10.^[-3:2]);
filename = sprintf('system_of_diffusion_eqns_1d_comp_time.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(4); clf;
plot(x_lo_res,u_FE_OTS_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
axis([-1 1 0.8 1.2]);
xlabel('x');
%title('Forward Euler OTS Solution')
filename = sprintf('system_of_diffusion_eqns_1d_u_FE_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(5); clf;
plot(x_lo_res,u_FE_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
axis([-1 1 0.8 1.2]);
xlabel('x');
%title('Forward Euler Solution')
filename = sprintf('system_of_diffusion_eqns_1d_u_FE_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(6); clf;
plot(x_lo_res,v_FE_OTS_lo_res,'bo')
hold on;
plot(x,v_exact,'r')
axis([-1 1 -1.5 1.5]);
xlabel('x');
%title('Forward Euler OTS Solution')
filename = sprintf('system_of_diffusion_eqns_1d_v_FE_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(7); clf;
plot(x_lo_res,v_FE_lo_res,'bo')
hold on;
plot(x,v_exact,'r')
axis([-1 1 -1.5 1.5]);
xlabel('x');
%title('Forward Euler Solution')
filename = sprintf('system_of_diffusion_eqns_1d_v_FE_soln.%s', print_suffix);
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

