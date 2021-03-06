%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for the 
% 1d Burgers equation with a single-hump solution for presentations.
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
print_format = 'png';
fig_dir = 'figures';
if ~exist(fig_dir, 'dir')
  mkdir(fig_dir);
end

% set flag for loading data from saved files (instead of recomputing solution)
use_saved_data = 1;
data_dir = 'data-burgers_eqn_1d-single_hump';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% simulation paramters
debug_on = 0;
timing_on = 1;

% physical parameters
nu = 0.1;  % viscosity
U  = 1.0;   % wave speed
R  = 10.0;  % effective Reynolds number

% time integration parameters
% NOTE: t_init = 0.0
t_final = 2.0;

% grid sizes to collect data on
grid_sizes = [50 100 200 400 800 1600 3200 6400];

% allocate memory for errors
err_FE_OTS = zeros(size(grid_sizes));
err_FE     = zeros(size(grid_sizes));

% allocate memory for computation times
comp_time_FE_OTS      = zeros(size(grid_sizes));
comp_time_FE          = zeros(size(grid_sizes));

% start clock for timing plot generation time
tic;


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
    dx = 10.0/N;
    dt_FE = dx^2/4/nu;

    % solve viscous Burgers equation using forward Euler with OTS
    disp('Forward Euler OTS');
    [u_FE_OTS, u_exact, x, timing_data_FE_OTS] = ...
      solveBurgersEqn1dForwardEulerOTS(nu, U, R, ...
                                       dx, ...
                                       t_final, ...
                                       debug_on, ...
                                       timing_on);

    % solve viscous Burgers equation using forward Euler
    disp('Forward Euler');
    [u_FE, u_exact, x, timing_data_FE] = ...
      solveBurgersEqn1dForwardEuler(nu, U, R, ...
                                    dx, dt_FE, ...
                                    t_final, ...
                                    debug_on, ...
                                    timing_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ...
         'u_FE_OTS', 'timing_data_FE_OTS', ...
         'u_FE', 'timing_data_FE', ...
         'u_exact', 'x');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  err = u_FE_OTS-u_exact;
  err_FE_OTS(i) = norm(err,'inf');
  err = u_FE-u_exact;
  err_FE(i) = norm(err,'inf');

% collect timing data
  comp_time_FE_OTS(i) = timing_data_FE_OTS(1);
  comp_time_FE(i) = timing_data_FE(1);

end % end loop over grid sizes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate solution at low grid resolution
% for illustration purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_lo_res = 100;
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
  dx = 10.0/N_lo_res;
  dt_FE = dx^2/4/nu;

  % solve viscous Burgers equation using forward Euler with OTS
  [u_FE_OTS_lo_res, u_exact_lo_res, x_lo_res] = ...
    solveBurgersEqn1dForwardEulerOTS(nu, U, R, ...
                                     dx, ...
                                     t_final, ...
                                     debug_on);

  % solve viscous Burgers equation using forward Euler
  [u_FE_lo_res, u_exact_lo_res, x_lo_res] = ...
    solveBurgersEqn1dForwardEuler(nu, U, R, ...
                                  dx, dt_FE, ...
                                  t_final, ...
                                  debug_on);

  % save solutions and timing data to MAT-files
  filename = 'data_lo_res';
  save([data_dir, '/', filename], ...
       'u_FE_OTS_lo_res', 'u_FE_lo_res', 'x_lo_res');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_FE_OTS = polyfit(log(grid_sizes),log(err_FE_OTS),1);
order_FE_OTS = -P_FE_OTS(1);
P_FE = polyfit(log(grid_sizes(2:end)),log(err_FE(2:end)),1);
order_FE = -P_FE(1);

P_comp_time_FE_OTS = polyfit(log(err_FE_OTS(5:end)), ...
                             log(comp_time_FE_OTS(5:end)),1);
P_comp_time_FE = polyfit(log(err_FE(5:end)), ...
                         log(comp_time_FE(5:end)),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [10 10000];
loglog(N_plot,exp(log(N_plot)*P_FE_OTS(1)+P_FE_OTS(2)),'k');
hold on;
plot(grid_sizes,err_FE_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS-NIDC)\nOrder = %1.1f', order_FE_OTS);
text(35,5e-8,order_str);

loglog(N_plot,exp(log(N_plot)*P_FE(1)+P_FE(2)),'k');
hold on;
plot(grid_sizes,err_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nOrder = %1.1f', order_FE);
text(400,1e-1,order_str);

axis([N_plot(1) N_plot(2) 1e-10 1e0]);
xlabel('N');
ylabel('L^\infty Error');
set(gca, 'YTick', 10.^[-10:2:0]);
set(gca, 'YMinorTick', 'off');
filename = sprintf('burgers_eqn_1d_error_vs_N.%s', print_format);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(2); clf;
err_plot = [1e-10 1];
loglog(err_plot,exp(log(err_plot)*P_comp_time_FE_OTS(1)+P_comp_time_FE_OTS(2)),'k');
hold on;
plot(err_FE_OTS, comp_time_FE_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS-NIDC)\nSlope = %1.2f', P_comp_time_FE_OTS(1));
text(1e-9,1e-2,order_str);

loglog(err_plot,exp(log(err_plot)*P_comp_time_FE(1)+P_comp_time_FE(2)),'k');
hold on;
plot(err_FE, comp_time_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nSlope = %1.1f', P_comp_time_FE(1));
text(4e-4,1e2,order_str);

axis([err_plot(1) err_plot(2) 1e-4 1e4]);
xlabel('L^\infty Error');
ylabel('Compute Time (s)');
set(gca, 'XTick', 10.^[-10:2:0]);
set(gca, 'YTick', 10.^[-4:2:4]);
filename = sprintf('burgers_eqn_1d_comp_time.%s', print_format);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(3); clf;
plot(x_lo_res,u_FE_OTS_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
axis([0 10 1 2]);
xlabel('x');
%title('Forward Euler OTS Solution')
filename = sprintf('burgers_eqn_1d_FE_OTS_soln.%s', print_format);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);

figure(4); clf;
plot(x_lo_res,u_FE_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
axis([0 10 1 2]);
xlabel('x');
%title('Forward Euler Solution')
filename = sprintf('burgers_eqn_1d_FE_soln.%s', print_format);
format_str = sprintf('-d%s',print_format);
print([fig_dir, '/', filename], format_str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display time for generating plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_total = toc;
disp('----------------------------------------');
disp_str = sprintf('Plot Generation Time: %f', t_total);
disp(disp_str);
disp('----------------------------------------');

