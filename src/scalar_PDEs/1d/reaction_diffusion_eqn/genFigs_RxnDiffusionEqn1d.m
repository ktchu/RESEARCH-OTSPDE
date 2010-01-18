%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for the 
% reaction-diffusion equation. 
%  
% Kevin T. Chu
% 2009 April
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
data_dir = 'data-rxn_diffusion_eqn_1d';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% set simulation parameters
debug_on = 0;
timing_on = 1;

% time integration parameters
t_final = 5.0;

% grid parameters
x_lo = -5;
x_hi =  5;

% grid sizes to collect data on
grid_sizes = [25 50 100 200 400 800];

% allocate memory for errors
err_FE_OTS         = zeros(size(grid_sizes));
err_FE             = zeros(size(grid_sizes));
err_NSFD           = zeros(size(grid_sizes));

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
    dt_FE = dx^2/4;

    % solve diffusion equation using forward Euler with OTS
    disp('Forward Euler OTS');
    [u_FE_OTS, u_exact, x, timing_data_FE_OTS] = ...
      solveRxnDiffusionEqn1dForwardEulerOTS(dx, t_final, ...
                                            debug_on, timing_on);

    % solve reaction-diffusion equation using forward Euler
    disp('Forward Euler');
    [u_FE, u_exact, x, timing_data_FE] = ...
      solveRxnDiffusionEqn1dForwardEuler(dx, dt_FE, t_final, ...
                                         debug_on, timing_on);

    % solve reaction-diffusion equation using NSFD
%    [u_NSFD, u_exact, x, timing_data_NSFD] = ...
%      solveRxnDiffusionEqn1dNSFD(dx, dt_NSFD, t_final, debug_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ... 
         'u_FE_OTS', 'u_FE', 'u_exact', 'x', ...
         'timing_data_FE_OTS', 'timing_data_FE');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  err = u_FE_OTS-u_exact;
  err_FE_OTS(i) = norm(err,'inf');
  err = u_FE-u_exact;
  err_FE(i) = norm(err,'inf');

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
  dx = (x_hi-x_lo)/N_lo_res;
  dt_FE = dx^2/4;

  % solve diffusion equation using forward Euler with OTS
  disp('Forward Euler OTS');
  [u_FE_OTS_lo_res, u_exact_lo_res, x_lo_res] = ...
    solveRxnDiffusionEqn1dForwardEulerOTS(dx, t_final, debug_on);

  % solve diffusion equation using forward Euler
  disp('Forward Euler');
  [u_FE_lo_res, u_exact_lo_res, x_lo_res] = ...
    solveRxnDiffusionEqn1dForwardEuler(dx, dt_FE, t_final, debug_on);

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
order_str = sprintf('Forward Euler OTS\nOrder = %1.1f', order_FE_OTS);
text(15,1e-6,order_str);

loglog(N_plot,exp(log(N_plot)*P_FE(1)+P_FE(2)),'k');
hold on;
plot(grid_sizes,err_FE, 'bs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler\nOrder = %1.1f', order_FE);
text(200,5e-2,order_str);

axis([10 1000 1e-10 1]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('rxn_diffusion_eqn_1d_error_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(2); clf;
err_plot = [1e-10 1e0];
loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE_OTS(1)+P_comp_time_FE_OTS(2)), ...
       'k');
hold on;
loglog(err_FE_OTS(2:end), comp_time_FE_OTS(2:end), 'bo', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','b');
order_str = sprintf('Forward Euler OTS\nSlope = %1.1f', ...
  comp_time_exp_FE_OTS);
text(3e-10,5e-1,order_str);

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE(1)+P_comp_time_FE(2)), ...
       'k');
hold on;
loglog(err_FE(2:end), comp_time_FE(2:end), 'rs', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nSlope = %1.1f', comp_time_exp_FE);
text(5e-5,80,order_str);

axis([1e-10 1e-2 1e-2 1e4]);
xlabel('L^\infty Error');
ylabel('Compute Time (s)');
filename = sprintf('rxn_diffusion_eqn_1d_comp_time.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(3); clf;
plot(x_lo_res,u_FE_OTS_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
xlabel('x');
%title('Forward Euler OTS Solution')
filename = sprintf('rxn_diffusion_eqn_1d_FE_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(4); clf;
plot(x_lo_res,u_FE_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
xlabel('x');
%title('Forward Euler Solution')
filename = sprintf('rxn_diffusion_eqn_1d_FE_soln.%s', print_suffix);
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

