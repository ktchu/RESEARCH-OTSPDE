%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solve4thOrderParabolicEqn1dForwardEuler(), 
% solve4thOrderParabolicEqn1dForwardEulerOTS(),
% solve4thOrderParabolicEqn1dCrankNicholson2ndOrder(), and
% solve4thOrderParabolicEqn1dCrankNicholson4thOrder().
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
t_final = 0.5;
t_final = 0.001;

% set dx and dt
N = 50;
dx = 1/N;
dt_FE = dx^4*3/40;  % maximum stable time step for forward Euler
dt_FE = dx^4/16;
dt_CN_2ndOrder = dx/4;
dt_CN_4thOrder = dx^2/4;

% solve diffusion equation using forward Euler with OTS
debug_on = 0;
[u_FE_OTS, u_exact, x] = solve4thOrderParabolicEqn1dForwardEulerOTS( ...
                             dx, ...
                             t_init, t_final, ...
                             debug_on);

% solve diffusion equation using forward Euler
[u_FE, u_exact, x] = solve4thOrderParabolicEqn1dForwardEuler( ...
                             dx, dt_FE, ...
                             t_init, t_final, ...
                             debug_on);

% solve diffusion equation using Crank-Nicholson with 2nd-order spatial
% discretization
[u_CN_2ndOrder, u_exact, x] = ...
  solve4thOrderParabolicEqn1dCrankNicholson2ndOrder(dx, dt_CN_2ndOrder, ...
                                                    t_init, t_final, ...
                                                    debug_on);

% solve diffusion equation using Crank-Nicholson with 4th-order spatial
% discretization
[u_CN_4thOrder, u_exact, x] = ...
  solve4thOrderParabolicEqn1dCrankNicholson4thOrder(dx, dt_CN_4thOrder, ...
                                                    t_init, t_final, ...
                                                    debug_on);

% compute error
err_FE_OTS = u_FE_OTS-u_exact;
err_L_inf_FE_OTS = norm(err_FE_OTS,'inf')
err_FE = u_FE-u_exact;
err_L_inf_FE = norm(err_FE,'inf')
err_CN_2ndOrder = u_CN_2ndOrder-u_exact;
err_L_inf_CN_2ndOrder = norm(err_CN_2ndOrder,'inf')
err_CN_4thOrder = u_CN_4thOrder-u_exact;
err_L_inf_CN_4thOrder = norm(err_CN_4thOrder,'inf')


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
plot(x,u_CN_2ndOrder,'bo')
hold on;
plot(x,u_exact,'r')
title('2nd-Order Crank-Nicholson Solution')

figure(6); clf;
plot(x,err_CN_2ndOrder);
title('Error in 2nd-Order Crank-Nicholson Solution')

figure(7); clf;
plot(x,u_CN_4thOrder,'bo')
hold on;
plot(x,u_exact,'r')
title('4th-Order Crank-Nicholson Solution')

figure(8); clf;
plot(x,err_CN_4thOrder);
title('Error in 4th-Order Crank-Nicholson Solution')

