%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of
% solveRxnDiffEqnForwardEulerOTS2d() and
% solveRxnDiffEqnForwardEuler2d().
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
g0   = 1;
R_in = 2;

% time integration parameters
% NOTE: t_init = 0.0
t_final = 1.0;
t_final = 0.1;

% set dx and dt
N = 100;
x_lo = -10;
x_hi = 10;
dx = (x_hi-x_lo)/N;
dt_FE = dx^2/3;  % suboptimal time step

% set simulation parameters
debug_on  = 1;
timing_on = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve 2d diffusion equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve diffusion equation using forward Euler with OTS
disp('----------------------------------------------------------');
disp('Solving RxnDiff Eqn Using Forward Euler with OTS');
disp('----------------------------------------------------------');
[u_FE_OTS, u_exact, X, Y, timing_data] = ...
   solveRxnDiffEqnForwardEulerOTS2d(g0, R_in, ...
                                    dx, ...
                                    t_final, ...
                                    debug_on, timing_on);

disp('----------------------------------------------------------');
disp('Solving RxnDiff Eqn Using Forward Euler without OTS');
disp('----------------------------------------------------------');
[u_FE, u_exact, X, Y, timing_data] = ...
   solveRxnDiffEqnForwardEuler2d(g0, R_in, ...
                                 dx, dt_FE, ...
                                 t_final, ...
                                 debug_on, timing_on);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_FE_OTS = u_FE_OTS-u_exact;
err_L_inf_FE_OTS = norm(err_FE_OTS,'inf')
err_FE = u_FE-u_exact;
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
[Y_plot,X_plot] = meshgrid(x_plot,y_plot);


figure(1); clf;
u_plot = interp2(reshape(X,N+1,N+1), ...
                 reshape(Y,N+1,N+1), ...
                 reshape(u_FE_OTS,N+1,N+1), ...
                 X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,u_plot);
title_string = sprintf('Forward Euler OTS Solution (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');

% plot current error
figure(2); clf;
err_plot = interp2(reshape(X,N+1,N+1), ...
                   reshape(Y,N+1,N+1), ...
                   reshape(err_FE_OTS,N+1,N+1), ...
                   X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,err_plot);
title_string = sprintf('Forward Euler OTS Error (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');


figure(3); clf;
u_plot = interp2(reshape(X,N+1,N+1), ...
                 reshape(Y,N+1,N+1), ...
                 reshape(u_FE,N+1,N+1), ...
                 X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,u_plot);
title_string = sprintf('Forward Euler Solution (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');

% plot current error
figure(4); clf;
err_plot = interp2(reshape(X,N+1,N+1), ...
                   reshape(Y,N+1,N+1), ...
                   reshape(err_FE,N+1,N+1), ...
                   X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,err_plot);
title_string = sprintf('Forward Euler Error (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');

