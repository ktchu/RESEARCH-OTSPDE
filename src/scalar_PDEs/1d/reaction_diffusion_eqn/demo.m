%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solveRxnDiffusionEqn1dForwardEuler(),
% solveRxnDiffusionEqn1dForwardEulerOTS(), and
% solveRxnDiffusionEqn1dNSFD().
%  
% Kevin T. Chu
% 2009 April
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% time integration parameters
t_final = 5.0;

% set dx and dt
x_lo = -5;
x_hi =  5;
N = 50;
dx = (x_hi-x_lo)/N;
dt_FE = 0.25*dx^2;
dt_NSFD = 0.5*dx^2;

% debug flag
debug_on = 0;

% solve diffusion equation using forward Euler with OTS
[u_FE_OTS, u_exact, x] = solveRxnDiffusionEqn1dForwardEulerOTS( ...
                             dx, t_final, debug_on);

% solve reaction-diffusion equation using forward Euler
[u_FE, u_exact, x] = solveRxnDiffusionEqn1dForwardEuler(dx, dt_FE, t_final, ...
                                                        debug_on);

% solve reaction-diffusion equation using non-standard finite difference scheme
%[u_NSFD, u_exact, x] = solveRxnDiffusionEqn1dNSFD(dx, dt_NSFD, t_final, ...
%                                                  debug_on);

% compute error
err_FE_OTS = u_FE_OTS-u_exact;
err_L_inf_FE_OTS = norm(err_FE_OTS,'inf')
err_FE = u_FE-u_exact;
err_L_inf_FE = norm(err_FE,'inf')
%err_NSFD = u_NSFD-u_exact;
%err_L_inf_NSFD = norm(err_NSFD,'inf')

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

%figure(5); clf;
%plot(x,u_NSFD,'bo')
%hold on;
%plot(x,u_exact,'r')
%title('NSFD Solution')
%
%figure(6); clf;
%plot(x,err_NSFD);
%title('Error in NSFD Solution')
