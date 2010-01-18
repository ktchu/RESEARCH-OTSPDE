%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script computes the solutions of the 4th-order parabolic
% equation in two spatial dimensions
%
%   u_t = -(u_xxxx + 2 u_xxyy + u_yyyy) + f
%
% where
%
%   f(x,y,t) = -3*sin(2*pi*x)*sin(5*pi*y)*sin(3*t)
%            + 841*pi^4*sin(2*pi*x)*sin(5*pi*y)*cos(3*t)
%            + 10*pi^2*sin(pi*x)*sin(3*pi*y)*exp(-10*pi^2*t)
%            + 100*pi^4*sin(pi*x)*sin(3*pi*y)*(1-exp(-10*pi^2*t)
%
% on the domain -1 < x,y < 1 with homogeneous Dirichlet boundary conditions 
% imposed at all boundaries.  The numerical solution is computed on a 
% node-centered grid using forward Euler time integration with a 
% fourth-order central difference approximation to the bilaplacian.  The 
% optimal time step for this problem is dt = 7/120*dx^4.
%
% NOTE:
% - For the purposes of illustration, the boundary conditions for the 
%   example are imposed by filling the ghostcells from the exact 
%   solution.
%
% Kevin T. Chu
% 2007 October
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% Set timing_run to 1 to disable error calculation and plotting
% during time integration loop.
timing_run = 0;

% time integration parameters
t_init  = 0.0;
t_final = 0.01;

% number of grid cells to use in computation
% NOTE: number of grid points in computational domain is (N_cell+1)
N_cell = 50;

% maximum number of grid cells to use for plotting
N_plot = 101;

% construct grid
N = N_cell+7;
num_gridpts = N*N;
dx = 2/N_cell; dy = dx;
x = -1-3*dx:dx:1+3*dx; y = x;   % include ghostcells
[Y,X] = meshgrid(x,y);
X = reshape(X,num_gridpts,1);
Y = reshape(Y,num_gridpts,1);

% set dt
dt = 7*dx^4/120;  % optimal time step
dt = dx^4/3;  % suboptimal time step


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for main computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up plotting grid
if (N < N_plot)
  N_plot = N;
end
dx_plot = 2/(N_plot-1);
x_plot = 0:dx_plot:1; y_plot = x_plot;
[Y_plot,X_plot] = meshgrid(x_plot,y_plot);

% start clock to measure time to prepare for computation
disp('Constructing Bilaplacian operators ...');
if (timing_run == 1)
  tic;
end

% construct 2nd-order Biaplacian operators (with boundary conditions) 
e = ones(N,1);
L_D = 1/dx^2*spdiags([e -2*e e], -1:1, N, N);
Lx = kron(speye(N),L_1D);
Ly = kron(L_1D,speye(N));

% construct 5pt-Bilaplacian operator (with boundary conditions) 
L_5pt = Lx + Ly; 

% construct 9pt-Laplacian operator (with boundary conditions) 
L_9pt = L_5pt + dx*dx/6*Lx*Ly;
L_9pt(1:N:N*N,:) = 0;         % Dirichlet BC x = 0
L_9pt(N:N:N*N,:) = 0;         % Dirichlet BC x = 1
L_9pt(1:N,:) = 0;             % Dirichlet BC y = 0
L_9pt((N-1)*N+1:N*N,:) = 0;   % Dirichlet BC y = 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve heat equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% measure Laplacian construction time and
% restart clock to measure computation time
if (timing_run == 1)
  time_laplacian_construction = toc;
  tic;
end

% initialize temperature
disp('Initializing termperature field ...');
if (source_term_type == 1)
  T = 1 + Y/3 + X/2 + X.*Y/4 ...
    - sin(pi*X).*sin(pi*Y) + sin(2*pi*X).*sin(3*pi*Y);
elseif (source_term_type == 2)
  T = 1 + Y/3 + X/2 + X.*Y/4 + sin(5*pi*X).*sin(pi*Y);
elseif (source_term_type == 3) 
  T = 1 + Y/3 + X/2 + X.*Y/4 ...
    + sin(2*pi*X).*sin(3*pi*Y) - sin(5*pi*X).*sin(pi*Y);
else
  error('Invalid source term type');
end


% forward Euler time integration
disp('Solving heat equation ...');
t = t_init;
time_step = 1;
figure(1); clf; % clear figure
while (t < t_final)
 
  if (timing_run ~= 1)
 
    % compute exact solution
    if (source_term_type == 1)
      T_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
              - sin(pi*X).*sin(pi*Y)*exp(-2*kappa*pi^2*t) ...
              + sin(2*pi*X).*sin(3*pi*Y)*exp(-13*kappa*pi^2*t);
    elseif (source_term_type == 2)
      T_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
              - sin(pi*X).*sin(pi*Y)*(1-exp(-2*kappa*pi^2*t)) ...
              + sin(2*pi*X).*sin(3*pi*Y)*(1-exp(-13*kappa*pi^2*t)) ...
              + sin(5*pi*X).*sin(pi*Y)*(2-exp(-26*kappa*pi^2*t));
    elseif (source_term_type == 3)
      T_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
              + sin(2*pi*X).*sin(3*pi*Y)*exp(-5*kappa*pi^2*t) ...
              - sin(5*pi*X).*sin(pi*Y)*exp(-2*kappa*pi^2*t);
    else
      error('Invalid source term type');
    end
    err = T-T_exact;
    err_L_inf = norm(err,'inf');

    % plot current solution
    figure(1); clf;
    T_plot = interp2(reshape(X,N,N), ...
                     reshape(Y,N,N), ...
                     reshape(T,N,N), ...
                     X_plot,Y_plot,'*cubic');
    surf(x_plot,y_plot,T_plot);
    title_string = sprintf('t = %f',t);
    title(title_string);
    xlabel('x'); ylabel('y');

    % plot current error
    figure(2); clf;
    err_plot = interp2(reshape(X,N,N), ...
                       reshape(Y,N,N), ...
                       reshape(err,N,N), ...
                       X_plot,Y_plot,'*cubic');
    surf(x_plot,y_plot,err_plot);
    title_string = sprintf('t = %f',t);
    title(title_string);
    xlabel('x'); ylabel('y');
    drawnow

    % display current status of solution
    disp('------------------------------------------------------------------');
    status_str = ...
      sprintf('Time = %f, Time step = %d, L infinity Error = %g', ...
              t, time_step, err_L_inf);
    disp(status_str); 

  end

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t;
  end

  % impose boundary conditions
  T(1:N) = 1+x/2;  % y = 0
  T(num_gridpts-N+1:num_gridpts) = 4/3+3/4*x;  % y = 1
  T(1:N:num_gridpts) = 1+y/3;  % x = 0
  T(N:N:num_gridpts) = 1.5+7/12*y;  % x = 1

  % compute source terms
  if (source_term_type == 1)
    f = 0;
    f_t = 0;
  elseif (source_term_type == 2)
    f = -2*kappa*pi^2*sin(pi*X).*sin(pi*Y) ...
      + 13*kappa*pi^2*sin(2*pi*X).*sin(3*pi*Y) ...
      + 52*kappa*pi^2*sin(5*pi*X).*sin(pi*Y);
    f_t = 0;
  elseif (source_term_type == 3)
    f = 8*kappa*pi^2*sin(2*pi*X).*sin(3*pi*Y)*exp(-5*kappa*pi^2*t) ...
      - 24*kappa*pi^2*sin(5*pi*X).*sin(pi*Y)*exp(-2*kappa*pi^2*t);
    f_t = -40*kappa^2*pi^4*sin(2*pi*X).*sin(3*pi*Y)*exp(-5*kappa*pi^2*t) ...
        + 48*kappa^2*pi^4*sin(5*pi*X).*sin(pi*Y)*exp(-2*kappa*pi^2*t);
  else
    error('Invalid source term type');
  end

  % update solution
  if (source_term_type == 1)
    T = T + dt*kappa*(L_9pt*T);
  elseif (source_term_type == 2)
    T = T + dt*kappa*(L_9pt*T) + dt*f + 0.5*dt^2*kappa*L_5pt*f; 
  elseif (source_term_type == 3)
    T = T + dt*kappa*(L_9pt*T) + dt*f + 0.5*dt^2*(f_t + kappa*L_5pt*f); 
  else 
    error('Invalid source term selection');
  end

  % update time and time step
  t = t + dt;
  time_step = time_step + 1;

end

% measure time to solve heat equation
if (timing_run == 1)
  time_solve = toc;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute exact solution and error
if (source_term_type == 1)
  T_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
          - sin(pi*X).*sin(pi*Y)*exp(-2*kappa*pi^2*t) ...
          + sin(2*pi*X).*sin(3*pi*Y)*exp(-13*kappa*pi^2*t);
elseif (source_term_type == 2)
  T_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
          - sin(pi*X).*sin(pi*Y)*(1-exp(-2*kappa*pi^2*t)) ...
          + sin(2*pi*X).*sin(3*pi*Y)*(1-exp(-13*kappa*pi^2*t)) ...
          + sin(5*pi*X).*sin(pi*Y)*(2-exp(-26*kappa*pi^2*t));
elseif (source_term_type == 3)
  T_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
          + sin(2*pi*X).*sin(3*pi*Y)*exp(-5*kappa*pi^2*t) ...
          - sin(5*pi*X).*sin(pi*Y)*exp(-2*kappa*pi^2*t);
else
  error('Invalid source term type');
  error('Invalid source term type');
end
err = T-T_exact;
err_L_inf = norm(err,'inf');

% output results statistics
disp(' ');
disp('==================================================================');
disp('Computation Statistics');
if (timing_run == 1)
  lapl_construct_time_str = ...
    sprintf('  Laplacian Construction Time: %f', time_laplacian_construction);
  comp_time_str = sprintf('  Solution Time: %f', time_solve);
  disp(lapl_construct_time_str);
  disp(comp_time_str);
end
err_L_inf_str = sprintf('  L infinity Norm of Error: %g', err_L_inf);
disp(err_L_inf_str);
disp('==================================================================');

% impose boundary conditions
T(1:N) = 1+x/2;  % y = 0
T(num_gridpts-N+1:num_gridpts) = 4/3+3/4*x;  % y = 1
T(1:N:num_gridpts) = 1+y/3;  % x = 0
T(N:N:num_gridpts) = 1.5+7/12*y;  % x = 1

% plot final result
figure(1); clf;
T_plot = interp2(reshape(X,N,N), ...
                 reshape(Y,N,N), ...
                 reshape(T,N,N), ...
                 X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,T_plot);
title_string = sprintf('t = %f',t);
title(title_string);
xlabel('x'); ylabel('y');

% plot current error
figure(2); clf;
err_plot = interp2(reshape(X,N,N), ...
                   reshape(Y,N,N), ...
                   reshape(err,N,N), ...
                   X_plot,Y_plot,'*cubic');
surf(x_plot,y_plot,err_plot);
title_string = sprintf('t = %f',t);
title(title_string);
xlabel('x'); ylabel('y');
drawnow

