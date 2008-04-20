%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of
% solveAdvectionEqnForwardEulerOTS2d() and
% solveAdvectionEqnForwardEuler2d()
%
% Kevin T. Chu
% 2008 February
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% physical parameters
A_x = -1.0;  % flow speed in x-direction
A_y = -3.0;  % flow speed in y-direction

% time integration parameters
% NOTE: t_init = 0.0
t_final = 1.0;

% set dx and dt
N = 100;
x_lo = -10;
x_hi =  10;
dx = (x_hi-x_lo)/N;
dt_FE = 0.5*dx/( abs(A_x) + abs(A_y) );  % suboptimal time step

% set simulation parameters
debug_on  = 1;
timing_on = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve 2d advection equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve advection equation using forward Euler with OTS
disp('----------------------------------------------------------');
disp('Solving Advection Eqn Using Forward Euler with OTS');
disp('----------------------------------------------------------');
[u_FE_OTS, u_exact_FE_OTS, dy_FE_OTS, X_FE_OTS, Y_FE_OTS, timing_data] = ...
   solveAdvectionEqnForwardEulerOTS2d(A_x, A_y, ...
                                      dx, ...
                                      t_final, ...
                                      debug_on, timing_on);

disp('----------------------------------------------------------');
disp('Solving Advection Eqn Using Forward Euler without OTS');
disp('----------------------------------------------------------');
[u_FE, u_exact_FE, X_FE, Y_FE, timing_data] = ...
   solveAdvectionEqnForwardEuler2d(A_x, A_y, ...
                                   dx, dt_FE, ...
                                   t_final, ...
                                   debug_on, timing_on);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_FE_OTS = u_FE_OTS-u_exact_FE_OTS;
err_L_inf_FE_OTS = norm(err_FE_OTS,'inf')
err_FE = u_FE-u_exact_FE;
err_L_inf_FE = norm(err_FE,'inf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate plot grid
N_plot = 101;
if (N < N_plot)
  N_plot = N;
end

dx_plot = (x_hi-x_lo)/(N_plot-1);
x_plot = x_lo:dx_plot:x_hi; y_plot = x_plot;
[X_plot,Y_plot] = meshgrid(x_plot,y_plot);


figure(1); clf;
Ny_FE_OTS = floor(20/dy_FE_OTS) + 1;
u_plot = interp2(reshape(X_FE_OTS,Ny_FE_OTS,N+1), ...
                 reshape(Y_FE_OTS,Ny_FE_OTS,N+1), ...
                 reshape(u_FE_OTS,Ny_FE_OTS,N+1), ...
                 X_plot, Y_plot, '*nearest');
mesh(x_plot,y_plot,u_plot);
title_string = sprintf('Forward Euler OTS Solution (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');

% plot current error
figure(2); clf;
err_plot = interp2(reshape(X_FE_OTS,Ny_FE_OTS,N+1), ...
                 reshape(Y_FE_OTS,Ny_FE_OTS,N+1), ...
                 reshape(err_FE_OTS,Ny_FE_OTS,N+1), ...
                 X_plot, Y_plot, '*nearest');
mesh(x_plot,y_plot,err_plot);
title_string = sprintf('Forward Euler OTS Error (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');


figure(3); clf;
u_plot = interp2(reshape(X_FE,N+1,N+1), ...
                 reshape(Y_FE,N+1,N+1), ...
                 reshape(u_FE,N+1,N+1), ...
                 X_plot,Y_plot,'*cubic');
mesh(x_plot,y_plot,u_plot);
title_string = sprintf('Forward Euler Solution (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');

% plot current error
figure(4); clf;
err_plot = interp2(reshape(X_FE,N+1,N+1), ...
                   reshape(Y_FE,N+1,N+1), ...
                   reshape(err_FE,N+1,N+1), ...
                   X_plot,Y_plot,'*cubic');
mesh(x_plot,y_plot,err_plot);
title_string = sprintf('Forward Euler Error (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');


