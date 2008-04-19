%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of
% solveDiffusionEqnForwardEulerOTS2d() and solveDiffusionEqnForwardEuler2d() 
% for irregular domains.
%
% Kevin T. Chu
% 2008 February
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% physical parameters
D = 2.0;  % diffusion coefficient

% order of accuracy for ghostcells
bc_interpolant_order = 4;

% time integration parameters
t_init  = 0.0;
t_final = 0.01;
t_final = 0.001;

% set dt and dx
N = 100;
dx = 2.0/N;
dt_FE = dx^2/D/4;
dt_FE = dx^2/D/6;

% compute level set function that defines the domain
% NOTE:  phi is set to be a signed distance function to make
%        it easier to identify boundary points based on the 
%        value of phi
N_grid = N+1;
num_gridpts = N_grid*N_grid;
dy = dx;
x = -1:dx:1; y = x;
[X,Y] = meshgrid(x,y);   % grid created with y as fastest direction
X = reshape(X,num_gridpts,1);
Y = reshape(Y,num_gridpts,1);

disp('Generating level set function for boundary ...');
theta = atan2(Y,X);  
phi = sqrt(X.^2+Y.^2) - 0.8;                        % circle
phi = sqrt(X.^2+Y.^2) - 0.6*(1+0.5*cos(5*theta));   % starfish
phi = sqrt(X.^2+Y.^2) - 0.6*(1+0.5*sin(5*theta));   % starfish
phi = reshape(phi,N_grid,N_grid);
phi = computeDistanceFunction2d(phi,[dx,dy]);
phi = reshape(phi,num_gridpts,1);

% source term
source_term_type = 3;

% set simulation parameters
debug_on  = 1;
timing_on = 1;

% set zero_level_set_tol and extrap_tol
zero_level_set_tol = -1;
extrap_tol = -1;
extrap_tol = dx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve 2d diffusion equation on irreg. domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve diffusion equation using forward Euler with OTS
disp('----------------------------------------------------------');
disp('Solving Diffusion Eqn Using Forward Euler with OTS');
disp('----------------------------------------------------------');
[u_FE_OTS, u_exact, X, Y, timing_data, idx_ghostcells] = ...
   solveDiffusionEqnForwardEulerOTS2d(D, ...
                                      source_term_type, ...
                                      phi, ...
                                      dx, ...
                                      t_final, ...
                                      bc_interpolant_order, ...
                                      debug_on, timing_on, ...
                                      zero_level_set_tol, extrap_tol);

% solve diffusion equation using forward Euler witout OTS
disp('----------------------------------------------------------');
disp('Solving Diffusion Eqn Using Forward Euler');
disp('----------------------------------------------------------');
[u_FE, u_exact, X, Y, timing_data, idx_ghostcells] = ...
   solveDiffusionEqnForwardEuler2d(D, ...
                                   source_term_type, ...
                                   phi, ...
                                   dx, dt_FE, ...
                                   t_final, ...
                                   bc_interpolant_order, ...
                                   debug_on, timing_on, ...
                                   zero_level_set_tol, extrap_tol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_FE_OTS = u_FE_OTS-u_exact;
err_L_inf_FE_OTS = norm(err_FE_OTS,'inf')
err_FE = u_FE-u_exact;
err_L_inf_FE = norm(err_FE,'inf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate plot grid
N_plot = 101;
if (N < N_plot)
  N_plot = N;
end
dx_plot = 2/(N_plot-1);
x_plot = -1:dx_plot:1; y_plot = x_plot;
[X_plot,Y_plot] = meshgrid(x_plot,y_plot);


figure(1); clf;
u_plot = interp2(reshape(X,N_grid,N_grid), ...
                 reshape(Y,N_grid,N_grid), ...
                 reshape(u_FE_OTS,N_grid,N_grid), ...
                 X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,u_plot);
title_string = sprintf('Forward Euler OTS Solution (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');

% plot current error
figure(2); clf;
err_plot = interp2(reshape(X,N_grid,N_grid), ...
                   reshape(Y,N_grid,N_grid), ...
                   reshape(err_FE_OTS,N_grid,N_grid), ...
                   X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,abs(err_plot));
title_string = sprintf('Forward Euler OTS Error (t = %f)',t_final);
title(title_string);
xlabel('x'); ylabel('y');
view([0 90]);
hold on;
idx_bdry = idx_ghostcells{1};
idx_ghostcells_edge = idx_ghostcells{2};
idx_ghostcells_corner = idx_ghostcells{3};
plot(X(idx_ghostcells_edge),Y(idx_ghostcells_edge),'k+');
plot(X(idx_ghostcells_corner),Y(idx_ghostcells_corner),'mx');
plot(X(idx_bdry),Y(idx_bdry),'co');
axis square

% plot contours of level set function and location of ghostcells
idx_bdry = idx_ghostcells{1};
idx_ghostcells_edge = idx_ghostcells{2};
idx_ghostcells_corner = idx_ghostcells{3};
figure(3); clf
contourf(reshape(X,N_grid,N_grid), ...
                 reshape(Y,N_grid,N_grid), ...
                 reshape(phi,N_grid,N_grid), ...
                 -1:0.1:1);
hold on;
contour(reshape(X,N_grid,N_grid), ...
        reshape(Y,N_grid,N_grid), ...
        reshape(phi,N_grid,N_grid),[0 0], ...
        'LineColor','b','LineWidth',2);
plot(X(idx_ghostcells_edge),Y(idx_ghostcells_edge),'k+');
plot(X(idx_ghostcells_corner),Y(idx_ghostcells_corner),'mx');
plot(X(idx_bdry),Y(idx_bdry),'co');
axis square


