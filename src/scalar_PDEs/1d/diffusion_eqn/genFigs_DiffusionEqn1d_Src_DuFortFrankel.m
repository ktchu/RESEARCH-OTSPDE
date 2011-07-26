%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script generates figures of the results for the diffusion 
% equation with a source term when solved using the DuFort-Frankel 
% scheme.
%  
% Kevin T. Chu
% 2009 May
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
data_dir = 'data-diffusion_eqn_1d_src_DF';
if ~exist(data_dir, 'dir')
  mkdir(data_dir);
end

% set simulation parameters
debug_on = 0;
timing_on = 1;

% source term type
source_term_type = 1;

% physical parameters
D = 0.5;  % diffusion coefficient 

% boundary conditions 
u_0 = 1;
dudx_1 = 0.5;

% time integration parameters
t_init  = 0.0;
t_final = 0.1;

% grid sizes to collect data on
grid_sizes = [25 50 100 200 400 800];

% allocate memory for errors
err_DF_OTS = zeros(size(grid_sizes));
err_DF     = zeros(size(grid_sizes));

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
    dx = 1/N;
    dt_DF = dx^2/3/D;

    % solve diffusion equation using DuFort-Frankel with OTS
    disp('DuFort-Frankel OTS');
    [u_DF_OTS, u_exact, x, timing_data_DF_OTS] = ...
       solveDiffusionEqn1dDuFortFrankelOTS(D, ...
                                           source_term_type, ...
                                           u_0, dudx_1, ...
                                           dx, ...
                                           t_init, t_final, ...
                                           debug_on, timing_on);

    % solve diffusion equation using DuFort-Frankel
    disp('DuFort-Frankel');
    [u_DF, u_exact, x, timing_data_DF] = ...
       solveDiffusionEqn1dDuFortFrankel(D, ...
                                        source_term_type, ...
                                        u_0, dudx_1, ...
                                        dx, dt_DF, ...
                                        t_init, t_final, ...
                                        debug_on, timing_on);

    % save solutions and timing data to MAT-files
    filename = sprintf('data_%d', N);
    save([data_dir, '/', filename], ... 
         'u_DF_OTS', 'u_DF', 'u_exact', 'x', ...
         'timing_data_DF_OTS', 'timing_data_DF');

  end % end case:  (use_saved_data ~= 1) ==> recompute solutions

  % compute error
  err = u_DF_OTS-u_exact;
  err_DF_OTS(i) = norm(err,'inf');
  err = u_DF-u_exact;
  err_DF(i) = norm(err,'inf');

  % collect timing data
  comp_time_DF_OTS(i) = timing_data_DF_OTS;
  comp_time_DF(i) = timing_data_DF;

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
  dx = 1/N_lo_res;
  dt_DF = dx^2/3/D;

  % solve diffusion equation using DuFort-Frankel with OTS
  disp('DuFort-Frankel OTS');
  [u_DF_OTS_lo_res, u_exact_lo_res, x_lo_res] = ...
     solveDiffusionEqn1dDuFortFrankelOTS(D, ...
                                         source_term_type, ...
                                         u_0, dudx_1, ...
                                         dx, ...
                                         t_init, t_final, ...
                                         debug_on);

  % solve diffusion equation using DuFort-Frankel
  disp('DuFort-Frankel');
  [u_DF_lo_res, u_exact_lo_res, x_lo_res] = ...
     solveDiffusionEqn1dDuFortFrankel(D, ...
                                      source_term_type, ...
                                      u_0, dudx_1, ...
                                      dx, dt_DF, ...
                                      t_init, t_final, ...
                                      debug_on);

  % save solutions and timing data to MAT-files
  filename = 'data_lo_res';
  save([data_dir, '/', filename], ...
       'u_DF_OTS_lo_res', 'u_DF_lo_res', 'x_lo_res');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_DF_OTS = polyfit(log(grid_sizes),log(err_DF_OTS),1);
order_DF_OTS = -P_DF_OTS(1);
P_DF = polyfit(log(grid_sizes),log(err_DF),1);
order_DF = -P_DF(1);

P_comp_time_DF_OTS = polyfit(log(err_DF_OTS(2:end)), ...
                             log(comp_time_DF_OTS(2:end)),1);
comp_time_exp_DF_OTS = P_comp_time_DF_OTS(1);
P_comp_time_DF = polyfit(log(err_DF(2:end)),log(comp_time_DF(2:end)),1);
comp_time_exp_DF = P_comp_time_DF(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf;
N_plot = [10 1000];
loglog(N_plot,exp(log(N_plot)*P_DF_OTS(1)+P_DF_OTS(2)),'k');
hold on;
plot(grid_sizes,err_DF_OTS, 'bo', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','b');
order_str = sprintf('DuFort-Frankel (OTS-NIDC)\nOrder = %1.1f', order_DF_OTS);
text(20,1e-9,order_str);

loglog(N_plot,exp(log(N_plot)*P_DF(1)+P_DF(2)),'k');
hold on;
plot(grid_sizes,err_DF, 'rs', ...
     'MarkerSize',14, ...
     'MarkerFaceColor','r');
order_str = sprintf('DuFort-Frankel\nOrder = %1.1f', order_DF);
text(140,5e-3,order_str);

axis([10 1000 1e-12 1e0]);
set(gca, 'ytick', [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0]);
xlabel('N');
ylabel('L^\infty Error');
filename = sprintf('diffusion_eqn_1d_src_DF_error_vs_N.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(2); clf;
err_plot = [1e-12 1e0];
loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_DF_OTS(1)+P_comp_time_DF_OTS(2)), ...
       'k');
hold on;
loglog(err_DF_OTS, comp_time_DF_OTS, 'bo', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','b');
order_str = sprintf('DuFort-Frankel (OTS-NIDC)\nSlope = %1.1f', comp_time_exp_DF_OTS);
text(9e-10,5e-3,order_str);

loglog(err_plot, ...
       exp(log(err_plot)*P_comp_time_DF(1)+P_comp_time_DF(2)), ...
       'k');
hold on;
loglog(err_DF, comp_time_DF, 'rs', ...
       'MarkerSize',14, ...
       'MarkerFaceColor','r');
order_str = sprintf('DuFort-Frankel\nSlope = %1.1f', comp_time_exp_DF);
text(7e-6,20,order_str);

axis([1e-10 1e-2 1e-4 1e2]);
set(gca, 'xtick', [1e-10, 1e-8, 1e-6, 1e-4, 1e-2]);
xlabel('L^\infty Error');
ylabel('Compute Time (s)');
filename = sprintf('diffusion_eqn_1d_src_DF_comp_time.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);


figure(3); clf;
plot(x_lo_res,u_DF_OTS_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
axis([0 1 1 1.6]);
xlabel('x');
%title('DuFort-Frankel OTS Solution')
filename = sprintf('diffusion_eqn_1d_src_DF_OTS_soln.%s', print_suffix);
format_str = sprintf('-d%s', print_format);
print([fig_dir, '/', filename], format_str);

figure(4); clf;
plot(x_lo_res,u_DF_lo_res,'bo')
hold on;
plot(x,u_exact,'r')
axis([0 1 1 1.6]);
xlabel('x');
%title('DuFort-Frankel Solution')
filename = sprintf('diffusion_eqn_1d_src_DF_soln.%s', print_suffix);
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

