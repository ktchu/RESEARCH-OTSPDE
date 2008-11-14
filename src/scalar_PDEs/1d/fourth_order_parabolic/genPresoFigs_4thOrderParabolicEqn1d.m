%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for the 
% fourth-order parabolic equation for presentations.
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
print_format = 'png';
fig_dir = 'figures';
if ~exist(fig_dir, 'dir')
  mkdir(fig_dir);
end

% set flag for loading data from saved files (instead of recomputing solution)
use_saved_data = 1;
data_dir = 'data-fourth_order_parabolic_1d';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% simulation parameters
debug_on  = 0;
timing_on = 1;

% time integration parameters
t_init = 0.0;
t_final = 0.003;

% grid sizes to collect data on
grid_sizes = [25 50 100 200 400 800 1600 3200 6400 12800 25600];

% allocate memory for errors
err_FE_OTS      = zeros(size(grid_sizes));
err_FE          = zeros(size(grid_sizes));
err_CN_2ndOrder = zeros(size(grid_sizes));
err_CN_4thOrder = zeros(size(grid_sizes));

% allocate memory for computation times
comp_time_FE_OTS      = zeros(size(grid_sizes));
comp_time_FE          = zeros(size(grid_sizes));
comp_time_CN_2ndOrder = zeros(size(grid_sizes));
comp_time_CN_4thOrder = zeros(size(grid_sizes));

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
    dx = 2.0/N;
    dt_FE = dx^4/16;
    dt_CN_2ndOrder = dx/4;
    dt_CN_4thOrder = dx^2/4;

    % only run forward Euler calculations for up to 400 points
    if N <= 400
      % solve fourth-order parabolic equation using forward Euler with OTS
      disp('Forward Euler OTS');
      [u_FE_OTS, u_exact, x, timing_data_FE_OTS] = ...
         solve4thOrderParabolicEqn1dForwardEulerOTS(dx, ...
                                                    t_init, t_final, ...
                                                    debug_on, timing_on);

      % solve fourth-order parabolic equation using forward Euler
      disp('Forward Euler');
      [u_FE, u_exact, x, timing_data_FE] = ...
         solve4thOrderParabolicEqn1dForwardEuler(dx, dt_FE, ...
                                                 t_init, t_final, ...
                                                 debug_on, timing_on);
    end

    % only run 4th-order Crank-Nicholson calculations for up to 3200 points
    if N <= 3200
      % solve fourth-order parabolic equation using 4th-order Crank-Nicholson
      disp('4th-Order Crank-Nicholson');
      [u_CN_4thOrder, u_exact, x, timing_data_CN_4thOrder] = ...
         solve4thOrderParabolicEqn1dCrankNicholson4thOrder(dx, ...
                                                           dt_CN_4thOrder, ...
                                                           t_init, t_final, ...
                                                           debug_on, timing_on);
    end

    % solve fourth-order parabolic equation using 2nd-order Crank-Nicholson
    disp('2nd-Order Crank-Nicholson');
    [u_CN_2ndOrder, u_exact, x, timing_data_CN_2ndOrder] = ...
       solve4thOrderParabolicEqn1dCrankNicholson2ndOrder(dx, dt_CN_2ndOrder, ...
                                                         t_init, t_final, ...
                                                         debug_on, timing_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    if N <= 400
      save([data_dir, '/', filename], ...
           'u_FE_OTS', 'u_FE', ...
           'u_CN_2ndOrder', 'u_CN_4thOrder', ...
           'u_exact', 'x', ...
           'timing_data_FE_OTS', 'timing_data_FE', ...
           'timing_data_CN_2ndOrder', 'timing_data_CN_4thOrder');
    elseif N <= 3200
      save([data_dir, '/', filename], ...
           'u_CN_2ndOrder', 'u_CN_4thOrder', ...
           'u_exact', 'x', ...
           'timing_data_CN_2ndOrder', 'timing_data_CN_4thOrder');
    else
      save([data_dir, '/', filename], ...
           'u_CN_2ndOrder', 'u_exact', 'x', ...
           'timing_data_CN_2ndOrder');
    end

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  if N <= 400
    err = u_FE_OTS-u_exact;
    err_FE_OTS(i) = norm(err,'inf');
    err = u_FE-u_exact;
    err_FE(i) = norm(err,'inf');
  end

  if N <= 3200
    err = u_CN_4thOrder-u_exact;
    err_CN_4thOrder(i) = norm(err,'inf');
  end

  err = u_CN_2ndOrder-u_exact;
  err_CN_2ndOrder(i) = norm(err,'inf');

  % collect timing data
  if N <= 400
    comp_time_FE_OTS(i) = timing_data_FE_OTS(2);
    comp_time_FE(i) = timing_data_FE(2);
  end
  if N <= 3200
    comp_time_CN_4thOrder(i) = timing_data_CN_4thOrder(2);
  end
  comp_time_CN_2ndOrder(i) = timing_data_CN_2ndOrder(2);

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
  dt_FE = dx^4/16;
  dt_CN_2ndOrder = dx/4;
  dt_CN_4thOrder = dx^2/4;

  % solve fourth-order parabolic equation using forward Euler with OTS
  disp('Forward Euler OTS');
  [u_FE_OTS_lo_res, u_exact_lo_res, x_lo_res] = ...
     solve4thOrderParabolicEqn1dForwardEulerOTS(dx, ...
                                                t_init, t_final, ...
                                                debug_on);

  % solve fourth-order parabolic equation using forward Euler
  disp('Forward Euler');
  [u_FE_lo_res, u_exact_lo_res, x_lo_res] = ...
     solve4thOrderParabolicEqn1dForwardEuler(dx, dt_FE, ...
                                             t_init, t_final, ...
                                             debug_on);

  % solve fourth-order parabolic equation using 2nd-order Crank-Nicholson
  disp('2nd-Order Crank-Nicholson');
  [u_CN_2ndOrder_lo_res, u_exact_lo_res, x_lo_res] = ...
     solve4thOrderParabolicEqn1dCrankNicholson2ndOrder(dx, dt_CN_2ndOrder, ...
                                                       t_init, t_final, ...
                                                       debug_on, timing_on);

  % solve fourth-order parabolic equation using 4th-order Crank-Nicholson
  disp('4th-Order Crank-Nicholson');
  [u_CN_4thOrder_lo_res, u_exact_lo_res, x_lo_res] = ...
     solve4thOrderParabolicEqn1dCrankNicholson4thOrder(dx, dt_CN_4thOrder, ...
                                                       t_init, t_final, ...
                                                       debug_on, timing_on);

  % save solutions and timing data to MAT-files
  filename = 'data_lo_res';
  save([data_dir, '/', filename], ...
       'u_FE_OTS_lo_res', 'u_FE_lo_res', ...
       'u_CN_2ndOrder_lo_res', 'u_CN_4thOrder_lo_res',  ...
       'x_lo_res');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_FE_OTS = polyfit(log(grid_sizes(1:5)),log(err_FE_OTS(1:5)),1);
order_FE_OTS = -P_FE_OTS(1);
P_FE = polyfit(log(grid_sizes(1:5)),log(err_FE(1:5)),1);
order_FE = -P_FE(1);
P_CN_4thOrder = polyfit(log(grid_sizes(1:6)),log(err_CN_4thOrder(1:6)),1);
order_CN_4thOrder = -P_CN_4thOrder(1);
P_CN_2ndOrder = polyfit(log(grid_sizes(1:9)),log(err_CN_2ndOrder(1:9)),1);
order_CN_2ndOrder = -P_CN_2ndOrder(1);

P_comp_time_FE_OTS = polyfit(log(err_FE_OTS(1:5)), ...
                             log(comp_time_FE_OTS(1:5)),1);
comp_time_exp_FE_OTS = P_comp_time_FE_OTS(1);
P_comp_time_FE = polyfit(log(err_FE(1:5)),log(comp_time_FE(1:5)),1);
comp_time_exp_FE = P_comp_time_FE(1);
P_comp_time_CN_4thOrder = polyfit(log(err_CN_4thOrder(4:5)), ...
                                  log(comp_time_CN_4thOrder(4:5)),1);
comp_time_exp_CN_4thOrder = P_comp_time_CN_4thOrder(1);
P_comp_time_CN_2ndOrder = polyfit(log(err_CN_2ndOrder(7:9)), ...
                                  log(comp_time_CN_2ndOrder(7:9)),1);
comp_time_exp_CN_2ndOrder = P_comp_time_CN_2ndOrder(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure(1); clf;
N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_FE_OTS(1)+P_FE_OTS(2)),'k');
hold on;
plot(grid_sizes,err_FE_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS)\nOrder = %1.1f', order_FE_OTS);
text(30,1e-9,order_str);
annotation(f, 'arrow', [.63 .7], [.2 .2], ...
           'linewidth', 4, 'headstyle', 'plain');

N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_FE(1)+P_FE(2)),'k');
hold on;
plot(grid_sizes,err_FE, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nOrder = %1.1f', order_FE);
text(12,8e-6,order_str);
annotation(f, 'arrow', [.25 .37], [.55 .7], ...
           'linewidth', 4, 'headstyle', 'plain');

N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_CN_4thOrder(1)+P_CN_4thOrder(2)),'k');
hold on;
plot(grid_sizes,err_CN_4thOrder, 'cd', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','c');
order_str = sprintf('4th-Order Crank-Nicholson \nOrder = %1.1f', ...
                    order_CN_4thOrder);
text(12,7e-8,order_str);
annotation(f, 'arrow', [.4 .58], [.42 .55], ...
           'linewidth', 4, 'headstyle', 'plain');

N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_CN_2ndOrder(1)+P_CN_2ndOrder(2)),'k');
hold on;
plot(grid_sizes,err_CN_2ndOrder, 'g^', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','g');
order_str = sprintf('2nd-Order Crank-Nicholson \nOrder = %1.1f', ...
                    order_CN_2ndOrder);
text(80,7e-1,order_str);
annotation(f, 'arrow', [.7 .7], [.83 .75], ...
           'linewidth', 4, 'headstyle', 'plain');

axis([10 1000 1e-10 1e1]);
set(gca, 'ytick', [1e-10 1e-8 1e-6 1e-4 1e-2 1e0]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('4th_order_parabolic_eqn_1d_error_vs_N.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(2); clf;
N_plot = [10 100000];
loglog(N_plot,exp(log(N_plot)*P_CN_4thOrder(1)+P_CN_4thOrder(2)),'k');
hold on;
plot(grid_sizes,err_CN_4thOrder, 'cd', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','c');
order_str = sprintf('4th-Order Crank-Nicholson');
text(15,1e-9,order_str);

N_plot = [10 100000];
loglog(N_plot,exp(log(N_plot)*P_CN_2ndOrder(1)+P_CN_2ndOrder(2)),'k');
hold on;
plot(grid_sizes,err_CN_2ndOrder, 'g^', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','g');
order_str = sprintf('2nd-Order Crank-Nicholson');
text(400,5e-2,order_str);

axis([10 100000 1e-10 1e0]);
set(gca, 'ytick', [1e-10 1e-8 1e-6 1e-4 1e-2 1e0]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('4th_order_parabolic_eqn_1d_error_vs_N_roundoff.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


f = figure(3); clf;
err_plot = [1e-13 1e0];
loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE_OTS(1)+P_comp_time_FE_OTS(2)), ...
       'k');
hold on;
loglog(err_FE_OTS, comp_time_FE_OTS, 'bo', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','b');
order_str = sprintf('Forward Euler (OTS)\nSlope = %1.2f', comp_time_exp_FE_OTS);
text(5e-5,5e4,order_str);
annotation(f, 'arrow', [.58 .32], [.86 .78], ...
           'linewidth', 4, 'headstyle', 'plain');

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_FE(1)+P_comp_time_FE(2)), ...
       'k');
hold on;
loglog(err_FE, comp_time_FE, 'rs', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','r');
order_str = sprintf('Forward Euler\nSlope = %1.1f', comp_time_exp_FE);
text(1e-3,5e2,order_str);
annotation(f, 'arrow', [.78 .78], [.63 .52], ...
           'linewidth', 4, 'headstyle', 'plain');

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_CN_4thOrder(1) ...
          +P_comp_time_CN_4thOrder(2)), 'k');
hold on;
loglog(err_CN_4thOrder, comp_time_CN_4thOrder, 'cd', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','c');
order_str = sprintf('4th-Order Crank-Nicholson\nSlope = %1.2f', ...
                    comp_time_exp_CN_4thOrder);
text(1e-9,3e-4,order_str);
annotation(f, 'arrow', [.57 .60], [.27 .32], ...
           'linewidth', 4, 'headstyle', 'plain');

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_CN_2ndOrder(1) ...
          +P_comp_time_CN_2ndOrder(2)), ...
       'k');
hold on;
loglog(err_CN_2ndOrder, comp_time_CN_2ndOrder, 'g^', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','g');
order_str = sprintf('2nd-Order Crank-Nicholson\nSlope = %1.1f', ...
                    comp_time_exp_CN_2ndOrder);
text(4e-11,9e-3,order_str);
annotation(f, 'arrow', [.3 .35], [.39 .61], ...
           'linewidth', 4, 'headstyle', 'plain');

axis([2e-11 1e0 5e-5 1e6]);
set(gca, 'xtick', [1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0]);
set(gca, 'ytick', [1e-4 1e-2 1e0 1e2 1e4 1e6]);
xlabel('L^\infty Error');
ylabel('Compute Time (s)');
filename = sprintf('4th_order_parabolic_eqn_1d_comp_time.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(4); clf;
plot(x_lo_res,u_FE_OTS_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
xlabel('x');
%title('Forward Euler OTS Solution')
filename = sprintf('4th_order_parabolic_eqn_1d_FE_OTS_soln.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(5); clf;
plot(x_lo_res,u_FE_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
xlabel('x');
%title('Forward Euler Solution')
filename = sprintf('4th_order_parabolic_eqn_1d_FE_soln.%s', print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(6); clf;
plot(x_lo_res,u_CN_2ndOrder_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
xlabel('x');
%title('2nd-Order Crank-Nicholson Solution')
filename = sprintf('4th_order_parabolic_eqn_1d_CN_2ndOrder_soln.%s', ...
                   print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(7); clf;
plot(x_lo_res,u_CN_4thOrder_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
xlabel('x');
%title('4th-Order Crank-Nichsolson Solution')
filename = sprintf('4th_order_parabolic_eqn_1d_CN_4thOrder_soln.%s', ...
                   print_format);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display time for generating plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_total = cputime - t_start;
disp('----------------------------------------');
disp_str = sprintf('Plot Generation Time: %f', t_total);
disp(disp_str);
disp('----------------------------------------');
