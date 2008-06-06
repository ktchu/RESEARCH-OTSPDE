%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of
% solveDiffusionEqnForwardEulerOTS2d(),
% solveDiffusionEqnForwardEuler2d(), and
% solveDiffusionEqnCrankNicholson2d().
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
D = 0.25;  % diffusion coefficient

% time integration parameters
% NOTE: t_init = 0.0
t_final = 0.5;

% set dx and dt
N = 50;
dx = 1.0/N;
dt_FE = dx^2/D/3;  % suboptimal time step
dt_FE = 3/8*dx^2/D;  % suboptimal time step (maximum stable time step)
dt_FE = (3/8+0.01)*dx^2/D;  % suboptimal time step
dt_CN = dx/16/D;

% source term type
source_term_type = 3;

% set simulation parameters
debug_on  = 0;
timing_on = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve 2d diffusion equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve diffusion equation using forward Euler with OTS
disp('----------------------------------------------------------');
disp('Solving Diffusion Eqn Using Forward Euler with OTS');
disp('----------------------------------------------------------');
[u_FE_OTS, u_exact, X, Y, timing_data] = ...
   solveDiffusionEqnForwardEulerOTS2d(D, ...
                                      source_term_type, ...
                                      dx, ...
                                      t_final, ...
                                      debug_on, timing_on);

disp('----------------------------------------------------------');
disp('Solving Diffusion Eqn Using Forward Euler without OTS');
disp('----------------------------------------------------------');
[u_FE, u_exact, X, Y, timing_data] = ...
   solveDiffusionEqnForwardEuler2d(D, ...
                                   source_term_type, ...
                                   dx, dt_FE, ...
                                   t_final, ...
                                   debug_on, timing_on);

disp('----------------------------------------------------------');
disp('Solving Diffusion Eqn Using Crank-Nicholson');
disp('----------------------------------------------------------');
[u_CN, u_exact, X, Y, timing_data] = ...
   solveDiffusionEqnCrankNicholson2d(D, ...
                                     source_term_type, ...
                                     dx, dt_CN, ...
                                     t_final, ...
                                     debug_on, timing_on);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_FE_OTS = u_FE_OTS-u_exact;
err_L_inf_FE_OTS = norm(err_FE_OTS,'inf')
err_FE = u_FE-u_exact;
err_L_inf_FE = norm(err_FE,'inf')
err_CN = u_CN-u_exact;
err_L_inf_CN = norm(err_CN,'inf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate plot grid
N_plot = 51;
if (N+1 < N_plot)
  N_plot = N+1;
end

dx_plot = 1/(N_plot-1);
x_plot = 0:dx_plot:1; y_plot = x_plot;
[X_plot,Y_plot] = meshgrid(x_plot,y_plot);


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


figure(5); clf;
u_plot = interp2(reshape(X,N+1,N+1), ...
                 reshape(Y,N+1,N+1), ...
                 reshape(u_CN,N+1,N+1), ...
                 X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,u_plot);
title_string = sprintf('Crank-Nicholson Solution (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');

% plot current error
figure(6); clf;
err_plot = interp2(reshape(X,N+1,N+1), ...
                   reshape(Y,N+1,N+1), ...
                   reshape(err_CN,N+1,N+1), ...
                   X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,err_plot);
title_string = sprintf('Crank-Nicholson Error (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');

