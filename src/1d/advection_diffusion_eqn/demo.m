%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solveAdvDiffEqnForwardEuler1d() and solveAdvDiffEqnForwardEulerOTS1d().
%  
% Kevin T. Chu
% 2008 February
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% physical parameters
D = 5.0;  % diffusion coefficient 
A = 5.0;  % flow speed

% time integration parameters
% NOTE: t_init = 0.0
t_final = 1.0;

% set dx and dt
N = 100;
dx = 20.0/N;
dt_FE = dx^2/D/4;

% solve diffusion equation using forward Euler with OTS
debug_on = 1;
[u_FE_OTS, u_exact, x] = solveAdvDiffEqnForwardEulerOTS1d( ...
                             D, A, ...
                             dx, ...
                             t_final, ...
                             debug_on);

% solve diffusion equation using forward Euler
[u_FE, u_exact, x] = solveAdvDiffEqnForwardEuler1d( ...
                             D, A, ...
                             dx, dt_FE, ...
                             t_final, ...
                             debug_on);

% compute error
err_FE_OTS = u_FE_OTS-u_exact;
err_L_inf_FE_OTS = norm(err_FE_OTS,'inf')
err_FE = u_FE-u_exact;
err_L_inf_FE = norm(err_FE,'inf')

% plot results
figure(1); clf;
plot(x,u_FE_OTS,'bo')
hold on;
plot(x,u_exact,'r')
title('Forward Euler OTS Solution')

figure(2); clf;
plot(x,err_FE_OTS);
title('Error in Forward Euler OTS Solution')

figure(3); clf;
plot(x,u_FE,'bo')
hold on;
plot(x,u_exact,'r')
title('Forward Euler Solution')

figure(4); clf;
plot(x,err_FE);
title('Error in Forward Euler Solution')
