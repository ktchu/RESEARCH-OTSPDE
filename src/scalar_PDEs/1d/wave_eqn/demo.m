%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solveWaveEqn1dKPY() and solveWaveEqn1dKPY_OTS()
%  
% Kevin T. Chu
% 2008 August
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% physical parameters
c = 1.0;  % wave speed
c = 0.5;  % wave speed
use_source_term = 1;

% time integration parameters
% NOTE: t_init = 0.0
t_final = 3.0;

% set dx and dt
x_lo = -0.5;
x_hi =  0.5;
N = 200;
dx = (x_hi-x_lo)/N;
dt_KPY = dx/c;
dt_KPY = 0.5*dx/c;

% solve diffusion equation using forward Euler with OTS
debug_on = 0;
[u_KPY_OTS, u_exact, x] = solveWaveEqn1dKPY_OTS(c, ...
                                                use_source_term, ...
                                                N, ...
                                                t_final, ...
                                                debug_on);

% solve diffusion equation using forward Euler
[u_KPY, u_exact, x] = solveWaveEqn1dKPY(c, ...
                                        use_source_term, ...
                                        N, dt_KPY, ...
                                        t_final, ...
                                        debug_on);

% compute error
err_KPY_OTS = u_KPY_OTS-u_exact;
err_L_inf_KPY_OTS = norm(err_KPY_OTS,'inf')
err_KPY = u_KPY-u_exact;
err_L_inf_KPY = norm(err_KPY,'inf')

% plot results
figure(1); clf;
plot(x,u_KPY_OTS,'bo')
hold on;
plot(x,u_exact,'r')
title('Forward Euler OTS Solution')

figure(2); clf;
plot(x,err_KPY_OTS);
title('Error in Forward Euler OTS Solution')

figure(3); clf;
plot(x,u_KPY,'bo')
hold on;
plot(x,u_exact,'r')
title('Forward Euler Solution')

figure(4); clf;
plot(x,err_KPY);
title('Error in Forward Euler Solution')