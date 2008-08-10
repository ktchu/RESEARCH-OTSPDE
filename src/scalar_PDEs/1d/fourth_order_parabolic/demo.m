%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solve4thOrderParabolicEqnForwardEuler1d() and
% solve4thOrderParabolicEqnForwardEulerOTS1d(). 
%  
% Kevin T. Chu
% 2008 February
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% time integration parameters
t_init  = 0.0;
t_final = 0.0001;
t_final = 1e-6;

% set dx and dt
N = 100;
dx = 1/N;
dt_FE = dx^4*3/40;  % maximum stable time step for forward Euler
dt_FE = dx^4/16;

% solve diffusion equation using forward Euler with OTS
debug_on = 0;
[u_FE_OTS, u_exact, x] = solve4thOrderParabolicEqnForwardEulerOTS1d( ...
                             dx, ...
                             t_init, t_final, ...
                             debug_on);

% solve diffusion equation using forward Euler
[u_FE, u_exact, x] = solve4thOrderParabolicEqnForwardEuler1d( ...
                             dx, dt_FE, ...
                             t_init, t_final, ...
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

