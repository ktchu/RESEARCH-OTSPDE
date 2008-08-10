%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveRxnDiffEqnCrankNicholson2d() computes the solutions of the 2d
% reaction-diffusion equation
%
%   u_t = (u_xx + u_yy) - u ln u
%
% on the domain -10 < x,y < 10 with initial conditions
%
%   u(x,y) = g0 * exp( -(x^2+y^2)/R_in^2 )
%
% and boundary conditions given by
%
%   u(x,y,t) = exp( -(r/R(t))^2 )
%            * exp( exp(-t) * (ln g0 - 2*sigma/(4+R_in^2) ln (R(t)/R_in)) ) 
%
% where R(t) = sqrt((4+R_in^2) exp(t) - 4).  sigma = 2+2*(d-1) is a constant
% that depends on the dimension d of the problem.  These boundary conditions 
% are obtained from the analytical solution of the reaction-diffusion equation 
% derived in Petrovskii and Shigesada (2001).  This particular solution 
% corresponds to the case v = 1 (i.e. 2D problem).  The numerical solution is 
% computed on a node-centered grid using forward Euler time integration with 
% the standard second-order 9pt discretization of the Laplacian.
%
% USAGE:
%  function [u, u_exact, X, Y, timing_data] = ...
%    solveRxnDiffEqnForwardEuler2d(g0, ...
%                                  R_in, ...
%                                  dx, dt, ...
%                                  t_final, ...
%                                  debug_on, timing_on)
%
% Arguments:
% - g0:         initial value for g(t) in Petrovskii and Shigesada
% - R_in:       initial "radius" of concentration
% - dx:         grid spacing
% - dt:         time step
% - t_final:    final time
% - debug_on:   flag indicating whether debugging information
%               should be displayed.  To turn on debugging,
%               set debug_on to 1.
%               (default = 0)
% - timing_on:  flag indicating whether timing information
%               should be collected.  To activate timing,
%               set timing_on to 1.  
%               (default = 0)
%
% Return values:
% - u:                   numerical solution
% - u_exact:             analytical solution
% - X, Y:                grid points
% - timing_data:         array of the following form containing timing data
%                        [laplacian_construction time, solution time].
%                        If timing is not activated, timing_data is 
%                        set to [-1, -1].
%
% NOTES:
% - The grid spacing is assumed to be the same in both the x and y directions.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGE LOG:
% -----------
% 2008/02:  Initial version of code.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [u, u_exact, X, Y, timing_data] = ...
  solveRxnDiffEqnForwardEuler2d(g0, ...
                                R_in, ...
                                dx, dt, ...
                                t_final, ...
                                debug_on, timing_on)

% check arguments
if (nargin < 5)
  error('solveRxnDiffEqnForwardEuler2d: missing arguments');
end
if (nargin < 6)
  debug_on = 0;
end
if (nargin < 7)
  timing_on = 0;
end

% construct grid
dy = dx;
x_lo = -10; x_hi = 10;
N = (x_hi-x_lo)/dx+1;
num_gridpts = N*N;
x = x_lo:dx:x_hi; y = x;
[Y,X] = meshgrid(y,x);   % grid created with x as fastest direction
X = reshape(X,num_gridpts,1);
Y = reshape(Y,num_gridpts,1);
idx_bdry = unique([1:N, num_gridpts-N+1:num_gridpts, ...
                   1:N:num_gridpts, N:N:num_gridpts]);

% set maximum number of grid cells to use for debug plotting
if (debug_on)
  N_plot = 101;
  if (N < N_plot)
    N_plot = N;
  end

  dx_plot = (x_hi-x_lo)/(N_plot-1);
  x_plot = x_lo:dx_plot:x_hi; y_plot = x_plot;
  [Y_plot,X_plot] = meshgrid(y_plot,x_plot);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for main computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start clock to measure time to prepare for computation
disp('Constructing Laplacian operators ...');
if (timing_on == 1)
  t_start = cputime;
end

% construct Laplacian operators (with boundary conditions) 
e = ones(N,1);
L_1D = 1/dx^2*spdiags([e -2*e e], -1:1, N, N);
Lx = kron(speye(N),L_1D);
Ly = kron(L_1D,speye(N));

% construct 5pt-Laplacian operator (with boundary conditions) 
L_5pt = Lx + Ly; 
L_5pt(1:N:N*N,:) = 0;         % Dirichlet BC x = x_lo
L_5pt(N:N:N*N,:) = 0;         % Dirichlet BC x = x_hi
L_5pt(1:N,:) = 0;             % Dirichlet BC y = y_lo
L_5pt((N-1)*N+1:N*N,:) = 0;   % Dirichlet BC y = y_hi

% construct 9pt-Laplacian operator (with boundary conditions) 
L_9pt = L_5pt + dx*dx/6*Lx*Ly;
L_9pt(1:N:N*N,:) = 0;         % Dirichlet BC x = x_lo
L_9pt(N:N:N*N,:) = 0;         % Dirichlet BC x = x_hi
L_9pt(1:N,:) = 0;             % Dirichlet BC y = y_lo
L_9pt((N-1)*N+1:N*N,:) = 0;   % Dirichlet BC y = y_hi

% cache frequently used values
two_sigma = 2*(2+2*(1));
rSq = X.^2 + Y.^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve reaction-diffusion equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% measure Laplacian construction time and
% restart clock to measure computation time
if (timing_on == 1)
  time_laplacian_construction = cputime - t_start;
  t_start = cputime;
end

% initialize concentration
disp('Initializing concentration field ...');
u = g0 * exp(-rSq/R_in^2);

% forward Euler time integration
disp('Solving reaction-diffusion equation ...');
t = 0.0;
time_step = 1;
while (t < t_final)

  if (debug_on == 1)
 
    % compute exact solution and error
    R = sqrt((4+R_in^2)*exp(t) - 4);
    u_exact = exp(-rSq/R^2) ...
            * exp( exp(-t)*(log(g0) - two_sigma/(4+R_in^2) * log(R/R_in)) ); 
    err = u-u_exact;
    err_L_inf = norm(err,'inf');

    % plot current solution
    figure(1); clf;
    u_plot = interp2(reshape(X,N,N), ...
                     reshape(Y,N,N), ...
                     reshape(u,N,N), ...
                     X_plot,Y_plot,'*cubic');
    surf(x_plot,y_plot,u_plot);
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

  % update solution
  u = u + dt*(L_9pt*u) - dt*u.*log(u);

  % update time and time step
  t = t + dt;
  time_step = time_step + 1;

  % impose boundary conditions
  R = sqrt((4+R_in^2)*exp(t) - 4);
  u(idx_bdry) = exp(-rSq(idx_bdry)/R^2) ...
              * exp( exp(-t)*(log(g0) - two_sigma/(4+R_in^2) * log(R/R_in)) ); 

end

% output extra line when in debug mode for aesthetic purposes
if (debug_on)
  disp('------------------------------------------------------------------');
end

% measure time to solve diffusion equation
if (timing_on == 1)
  time_solve = cputime - t_start;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output timing statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (timing_on) 
  timing_data = [time_laplacian_construction, time_solve];

  if (debug_on)
    disp(' ');
    disp('==================================================================');
    disp('Computation Statistics');
    lapl_construct_time_str = ...
      sprintf('  Laplacian Construction Time: %f', time_laplacian_construction);
    comp_time_str = sprintf('  Solution Time: %f', time_solve);
    disp(lapl_construct_time_str);
    disp(comp_time_str);
    disp('==================================================================');
  end

else

  timing_data = [-1, -1];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute exact solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = sqrt((4+R_in^2)*exp(t) - 4);
u_exact = exp(-rSq/R^2) ...
        * exp( exp(-t)*(log(g0) - two_sigma/(4+R_in^2) * log(R/R_in)) ); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impose boundary conditions on numerical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u(idx_bdry) = exp(-rSq(idx_bdry)/R^2) ...
            * exp( exp(-t)*(log(g0) - two_sigma/(4+R_in^2) * log(R/R_in)) ); 
