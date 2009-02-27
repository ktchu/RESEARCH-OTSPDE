%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of solveDiffusionEqn1dFE(), 
% solveDiffusionEqn1dFE_OTS_OptGrid().
%  
% Kevin T. Chu
% 2009 February
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% simulation parameters
debug_on = 0;

% time integration parameters
% NOTE: t_init = 0.0
t_final = 0.1;

% set dx and dt
x_lo = 0.0;
x_hi = 1.0;
N = 25;
dx = (x_hi-x_lo)/N;
dt_FE = 0.5*dx^2/9;

% solve 1d diffusion equation using FE time integration with OTS 
% and optimal grid
[u_FE_OTS_OptGrid, u_exact_OTS_OptGrid, x_OTS_OptGrid] = ...
  solveDiffusionEqn1dFE_OTS_OptGrid(N, t_final, debug_on);

% solve 1d diffusion equation using FE time integration without OTS
[u_FE, u_exact, x] = solveDiffusionEqn1dFE(N, dt_FE, t_final, debug_on);


% compute error
err_FE = u_FE-u_exact;
err_L_inf_FE = norm(err_FE,'inf')
err_FE_OTS_OptGrid = u_FE_OTS_OptGrid-u_exact_OTS_OptGrid;
err_L_inf_FE_OTS_OptGrid = norm(err_FE_OTS_OptGrid,'inf')

% plot results
figure(1); clf;
plot(x,u_FE,'bo')
hold on;
plot(x,u_exact,'r')
title('FE Solution')
axis([x_lo x_hi 0 4]);

figure(2); clf;
plot(x,err_FE);
title('Error in FE Solution')

figure(3); clf;
plot(x_OTS_OptGrid,u_FE_OTS_OptGrid,'bo')
hold on;
plot(x_OTS_OptGrid,u_exact_OTS_OptGrid,'r')
title('FE-OTS-OptGrid Solution')
axis([x_lo x_hi 0 4]);

figure(4); clf;
plot(x_OTS_OptGrid,err_FE_OTS_OptGrid);
title('Error in FE-OTS-OptGrid Solution')

