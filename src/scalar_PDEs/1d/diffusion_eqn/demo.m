%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solveDiffusionEqnForwardEuler1d(),
% solveDiffusionEqnForwardEulerOTS1d(), and
% solveDiffusionEqnCrankNicholson1d().
%  
% Kevin T. Chu
% 2008 February
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% physical parameters
D = 1.0;  % diffusion coefficient 

% boundary conditions 
u_0 = 1;
dudx_1 = 0.5;

% time integration parameters
t_init  = 0.0;
t_final = 0.1;

% set dx and dt
N = 100;
dx = 1/N;
dt_FE = dx^2/4/D;
dt_CN = dx/8/D;

% source term type
source_term_type = 0;

% solve diffusion equation using forward Euler with OTS
debug_on = 0;
[u_FE_OTS, u_exact, x] = solveDiffusionEqnForwardEulerOTS1d( ...
                             D, ...
                             source_term_type, ...
                             u_0, dudx_1, ...
                             dx, ...
                             t_init, t_final, ...
                             debug_on);

% solve diffusion equation using forward Euler
[u_FE, u_exact, x] = solveDiffusionEqnForwardEuler1d( ...
                             D, ...
                             source_term_type, ...
                             u_0, dudx_1, ...
                             dx, dt_FE, ...
                             t_init, t_final, ...
                             debug_on);

% solve diffusion equation using Crank-Nicholson
[u_CN, u_exact, x] = solveDiffusionEqnCrankNicholson1d( ...
                             D, ...
                             source_term_type, ...
                             u_0, dudx_1, ...
                             dx, dt_CN, ...
                             t_init, t_final, ...
                             debug_on);

% compute error
err_FE_OTS = u_FE_OTS-u_exact;
err_L_inf_FE_OTS = norm(err_FE_OTS,'inf')
err_FE = u_FE-u_exact;
err_L_inf_FE = norm(err_FE,'inf')
err_CN = u_CN-u_exact;
err_L_inf_CN = norm(err_CN,'inf')

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

figure(5); clf;
plot(x,u_CN,'bo')
hold on;
plot(x,u_exact,'r')
title('Crank-Nicholson Solution')

figure(6); clf;
plot(x,err_CN);
title('Error in Crank-Nicholson Solution')
