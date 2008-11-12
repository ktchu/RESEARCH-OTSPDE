%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solveBurgersEqn1dForwardEuler() and solveBurgersEqn1dForwardEulerOTS().
%  
% Kevin T. Chu
% 2008 February
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% physical parameters
nu = 0.1;  % viscosity
U  = 1.0;   % wave speed
R  = 10.0;  % effective Reynolds number

% time integration parameters
% NOTE: t_init = 0.0
t_final = 2.0;

% set dx and dt
N = 400;
N = 200;
N = 100;
dx = 10.0/N;
dt_FE = dx^2/nu/4;

% solve diffusion equation using forward Euler with OTS
debug_on = 1;
[u_FE_OTS, u_exact, x] = solveBurgersEqn1dForwardEulerOTS( ...
                             nu, U, R, ...
                             dx, ...
                             t_final, ...
                             debug_on);

% solve diffusion equation using forward Euler
[u_FE, u_exact, x] = solveBurgersEqn1dForwardEuler( ...
                             nu, U, R, ...
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
