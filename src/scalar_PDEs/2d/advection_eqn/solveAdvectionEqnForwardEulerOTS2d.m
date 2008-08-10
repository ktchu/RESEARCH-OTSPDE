%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveAdvectionEqnForwardEulerOTS2d() computes the solutions of the 2d 
% advection equation
%
%   u_t + A_x*u_x + A_y*u_y = 0
%
% on the domain -\infty < x,y < \infty with initial conditions
% 
%   u(x,y) = 1 if |x| + |y| <= 2
%          = 0 if |x| + |y| >  2
%
% using forward Euler time integration with a first-order upwind 
% spatial discretization.  For computational purposes, the domain
% is truncated to -10 < x,y < 10 and the boundary conditions are taken
% is truncated to -10 < x,y < 10 and the boundary conditions are taken
% from the exact solution.  The optimal time step of dx/? and the optimal 
% grid spacing ratio of dy/dx = ?? are used.
%
% USAGE:
%  function [u, u_exact, dy, X, Y, timing_data] = ...
%    solveAdvectionEqnForwardEulerOTS2d(A_x, A_y, ...
%                                       dx, ...
%                                       t_final, ...
%                                       debug_on, timing_on)
%
% Arguments:
% - A_x, A_y:            components of the flow velocity
% - dx:                  grid spacing in the x-direction
% - t_final:             final time
% - debug_on:            flag indicating whether debugging information
%                        should be displayed.  To turn on debugging,
%                        set debug_on to 1.
%                        (default = 0)
% - timing_on:           flag indicating whether timing information
%                        should be collected.  To activate timing,
%                        set timing_on to 1.
%                        (default = 0)
%
% Return values:
% - u:                   numerical solution
% - u_exact:             analytical solution
% - dy:                  optimal grid spacing in the y-direction
% - X, Y:                grid points
% - timing_data:         array of the following form containing timing data
%                        [gradient_construction time, solution time].
%                        If timing is not activated, timing_data is
%                        set to [-1, -1].
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGE LOG:
% -----------
% 2008/02:  Changed code to be a function.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [u, u_exact, dy, X, Y, timing_data] = ...
  solveAdvectionEqnForwardEulerOTS2d(A_x, A_y, ...
                                     dx, ...
                                     t_final, ...
                                     debug_on, timing_on)

% check arguments
if (nargin < 4)
  error('solveAdvectionEqnForwardEulerOTS2d: missing arguments');
end
if (nargin < 5)
  debug_on = 0;
end
if (nargin < 6)
  timing_on = 0;
end

% compute grid spacing in y-direction
dy = dx*abs(A_y/A_x);

% construct grid
x_lo = -10;
x_hi =  10;
Nx = (x_hi-x_lo)/dx+1;
x = x_lo:dx:x_hi; 
y = x_lo:dy:x_hi; 
Ny = length(y);
num_gridpts = Nx*Ny;
[X,Y] = meshgrid(x,y);   % grid created with y as fastest direction
X = reshape(X,num_gridpts,1);
Y = reshape(Y,num_gridpts,1);
idx_bdry = unique([1:Nx, num_gridpts-Nx+1:num_gridpts, ...
                   1:Nx:num_gridpts, Nx:Nx:num_gridpts]);

% compute optimal time step 
dt = dx/(abs(A_x));

% set maximum number of grid cells to use for debug plotting
if (debug_on)
  N_plot = 51;

  dx_plot = (x_hi-x_lo)/(N_plot-1);
  x_plot = x_lo:dx_plot:x_hi; y_plot = x_plot;
  [X_plot,Y_plot] = meshgrid(x_plot,y_plot);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for main computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start clock to measure time to prepare for computation
disp('Constructing upwind gradient operators ...');
if (timing_on == 1)
  t_start = cputime;
end

% construct upwind gradient operators (with boundary conditions) 
e = ones(Nx,1);
if (A_x > 0) 
  G_1D_x = 1/dx*spdiags([-e e], [-1,0], Nx, Nx);
else
  G_1D_x = 1/dx*spdiags([-e e], [0,1], Nx, Nx);
end
G_1D_x(1,:) = 0;    % no need to compute gradient for grid point at x = -10
G_1D_x(end,:) = 0;  % no need to compute gradient for grid point at x = 10
G_x = kron(G_1D_x,speye(Ny));

e = ones(Ny,1);
if (A_y > 0) 
  G_1D_y = 1/dy*spdiags([-e e], [-1,0], Ny, Ny);
else
  G_1D_y = 1/dy*spdiags([-e e], [0,1], Ny, Ny);
end
G_1D_y(1,:) = 0;    % no need to compute gradient for grid point at y = -10
G_1D_y(end,:) = 0;  % no need to compute gradient for grid point at y = 10
G_y = kron(speye(Nx),G_1D_y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve advection equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% measure gradient construction time and
% restart clock to measure computation time
if (timing_on == 1)
  time_gradient_construction = cputime - t_start;
  t_start = cputime;
end

% initialize concentration
disp('Initializing concentration field ...');
u = zeros(size(X));
idx = find( abs(X) + abs(Y) <= 2 );
u(idx) = 1;

% forward Euler time integration
disp('Solving advection equation ...');
t = 0.0;
time_step = 1;
while (t < t_final)
 
  if (debug_on == 1)
 
    % compute exact solution and error
    u_exact = zeros(size(X));
    idx = find( abs(X-A_x*t) + abs(Y-A_y*t) <= 2 );
    u_exact(idx) = 1;

    err = u-u_exact;
    err_L_inf = norm(err,'inf');

    % plot current solution
    figure(1); clf;
    u_plot = interp2(reshape(X,Ny,Nx), ...
                     reshape(Y,Ny,Nx), ...
                     reshape(u,Ny,Nx), ...
                     X_plot,Y_plot,'*nearest');
    mesh(x_plot,y_plot,u_plot);
    title_string = sprintf('t = %f',t);
    title(title_string);
    xlabel('x'); ylabel('y');

    % plot current error
    figure(2); clf;
    err_plot = interp2(reshape(X,Ny,Nx), ...
                       reshape(Y,Ny,Nx), ...
                       reshape(err,Ny,Nx), ...
                       X_plot,Y_plot,'*nearest');
    mesh(x_plot,y_plot,err_plot);
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
  u = u - dt*(A_x*G_x + A_y*G_y)*u + dt^2*A_x*A_y*G_x*G_y*u;

  % update time and time step
  t = t + dt;
  time_step = time_step + 1;

  % update boundary conditions
  u(idx_bdry) = 0;
  idx = find( abs(X(idx_bdry)-A_x*t) + abs(Y(idx_bdry)-A_y*t) <= 2);
  u(idx_bdry(idx)) = 1;

end

% output extra line when in debug mode for aesthetic purposes
if (debug_on)
  disp('------------------------------------------------------------------');
end

% measure time to solve advection equation
if (timing_on == 1)
  time_solve = cputime - t_start;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output timing statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (timing_on) 
  timing_data = [time_gradient_construction, time_solve];

  if (debug_on)
    disp(' ');
    disp('==================================================================');
    disp('Computation Statistics');
    lapl_construct_time_str = ...
      sprintf('  Gradient Construction Time: %f', time_gradient_construction);
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
u_exact = zeros(size(X));
idx = find( abs(X-A_x*t) + abs(Y-A_y*t) <= 2 );
u_exact(idx) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impose boundary conditions on numerical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u(idx_bdry) = 0;
idx = find( abs(X(idx_bdry)-A_x*t) + abs(Y(idx_bdry)-A_y*t) <= 2);
u(idx_bdry(idx)) = 1;

