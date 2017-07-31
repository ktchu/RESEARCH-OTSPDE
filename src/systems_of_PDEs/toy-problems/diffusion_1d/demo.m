%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solveDiffusionEqn1dForwardEuler() and solveDiffusionEqn1dForwardEulerOTS().
%  
% Kevin T. Chu
% 2008 July
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% time integration parameters
t_init  = 0.0;
t_final = 0.1;

% set dx and dt
N = 100;
dx = 2/N;
dt_FE = dx^2/3;

% simulation parameters
debug_on = 0;

% solve diffusion equations using forward Euler with OTS
[u_FE_OTS, v_FE_OTS, u_exact, v_exact, x] = ...
  solveDiffusionEqn1dForwardEulerOTS(dx, ...
                                     t_init, t_final, ...
                                     debug_on);

% solve diffusion equations using forward Euler
[u_FE, v_FE, u_exact, v_exact, x] = ...
  solveDiffusionEqn1dForwardEuler(dx, dt_FE, ...
                                  t_init, t_final, ...
                                  debug_on);

% compute error
err_u_FE_OTS = u_FE_OTS-u_exact;
err_u_L_inf_FE_OTS = norm(err_u_FE_OTS,'inf')
err_v_FE_OTS = v_FE_OTS-v_exact;
err_v_L_inf_FE_OTS = norm(err_v_FE_OTS,'inf')
err_u_FE = u_FE-u_exact;
err_u_L_inf_FE = norm(err_u_FE,'inf')
err_v_FE = v_FE-v_exact;
err_v_L_inf_FE = norm(err_v_FE,'inf')

% plot results
figure(1); clf;
plot(x,u_FE_OTS,'bo')
hold on;
plot(x,u_exact,'r')
title('Forward Euler OTS Solution for u')

figure(2); clf;
plot(x,v_FE_OTS,'bo')
hold on;
plot(x,v_exact,'r')
title('Forward Euler OTS Solution for v')

figure(3); clf;
plot(x,err_u_FE_OTS);
title('Error in Forward Euler OTS Solution for u')

figure(4); clf;
plot(x,err_v_FE_OTS);
title('Error in Forward Euler OTS Solution for v')

figure(5); clf;
plot(x,u_FE,'bo')
hold on;
plot(x,u_exact,'r')
title('Forward Euler Solution for u')

figure(6); clf;
plot(x,v_FE,'bo')
hold on;
plot(x,v_exact,'r')
title('Forward Euler Solution for v')

figure(7); clf;
plot(x,err_u_FE);
title('Error in Forward Euler Solution for u')

figure(8); clf;
plot(x,err_v_FE);
title('Error in Forward Euler Solution for v')

