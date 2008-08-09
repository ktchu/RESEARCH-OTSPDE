%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solveRxnDiffusionEqn1dForwardEuler(),
% solveRxnDiffusionEqn1dForwardEulerNoInterp() and 
% solveRxnDiffusionEqn1dForwardEulerOTS().
%  
% Kevin T. Chu
% 2008 August
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% time integration parameters
t_init  = 0.0;
t_final = 1.0;

% physical parameters
D = 2;

% set dx and dt
N = 400;
x_lo = -10;
x_hi =  10;
dx = (x_hi-x_lo)/N;
dt_FE = dx^2/4/D;

% simulation parameters
debug_on = 0;

% solve reaction-diffusion equations using forward Euler with OTS
[u_FE_OTS, v_FE_OTS, u_exact, v_exact, x] = ...
  solveRxnDiffusionEqn1dForwardEulerOTS(D, ...
                                        N, ...
                                        t_init, t_final, ...
                                        debug_on);

% solve reaction-diffusion equations using forward Euler with correction
% terms and minimum OTS among equations but no interpolation in time.
[u_FE_NoInterp, v_FE_NoInterp, u_exact, v_exact, x] = ...
  solveRxnDiffusionEqn1dForwardEulerNoInterp(D, ...
                                             N, ...
                                             t_init, t_final, ...
                                             debug_on);

% solve reaction-diffusion equations using forward Euler
[u_FE, v_FE, u_exact, v_exact, x] = ...
  solveRxnDiffusionEqn1dForwardEuler(D, ...
                                     N, dt_FE, ...
                                     t_init, t_final, ...
                                     debug_on);

% compute error
err_u_FE_OTS = u_FE_OTS-u_exact;
err_u_L_inf_FE_OTS = norm(err_u_FE_OTS,'inf')
err_v_FE_OTS = v_FE_OTS-v_exact;
err_v_L_inf_FE_OTS = norm(err_v_FE_OTS,'inf')
err_u_FE_NoInterp = u_FE_NoInterp-u_exact;
err_u_L_inf_FE_NoInterp = norm(err_u_FE_NoInterp,'inf')
err_v_FE_NoInterp = v_FE_NoInterp-v_exact;
rr_v_L_inf_FE_NoInterp = norm(err_v_FE_NoInterp,'inf')
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
plot(x,u_FE_NoInterp,'bo')
hold on;
plot(x,u_exact,'r')
set(gca, 'XLim', [x_lo x_hi]);
title('Forward Euler (No Interp) Solution for u')

figure(6); clf;
plot(x,v_FE_NoInterp,'bo')
hold on;
plot(x,v_exact,'r')
set(gca, 'XLim', [x_lo x_hi]);
title('Forward Euler (No Interp) Solution for v')

figure(7); clf;
plot(x,err_u_FE_NoInterp);
set(gca, 'XLim', [x_lo x_hi]);
title('Error in Forward Euler (No Interp) Solution for u')

figure(8); clf;
plot(x,err_v_FE_NoInterp);
set(gca, 'XLim', [x_lo x_hi]);
title('Error in Forward Euler (No Interp) Solution for v')

figure(9); clf;
plot(x,u_FE,'bo')
hold on;
plot(x,u_exact,'r')
set(gca, 'XLim', [x_lo x_hi]);
title('Forward Euler Solution for u')

figure(10); clf;
plot(x,v_FE,'bo')
hold on;
plot(x,v_exact,'r')
set(gca, 'XLim', [x_lo x_hi]);
title('Forward Euler Solution for v')

figure(11); clf;
plot(x,err_u_FE);
set(gca, 'XLim', [x_lo x_hi]);
title('Error in Forward Euler Solution for u')

figure(12); clf;
plot(x,err_v_FE);
set(gca, 'XLim', [x_lo x_hi]);
title('Error in Forward Euler Solution for v')

