%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqnCrankNicholson2d() computes the solutions of the 2d
% diffusion equation
%
%   u_t = D (u_xx + u_yy) + f
%
% on the domain 0 < x,y < 1 with Dirichlet boundary conditions imposed 
% at all boundaries.  The numerical solution is computed on a 
% node-centered grid using Crank-Nicholson time integration using 
% a standard second-order 5pt discretization of the Laplacian.  
%
% USAGE:
%  function [u, u_exact, X, Y, timing_data] = ...
%    solveDiffusionEqnCrankNicholson2d(D, ...
%                                      source_term_type, ...
%                                      dx, dt, ...
%                                      t_final, ...
%                                      debug_on, timing_on)
%
% Arguments:
% - D:                   diffusion coefficient
% - source_term_type:    type of source term to use in computation.
%                        1:  f = 0
%
%                        2:  f = -2*D*pi^2*sin(pi*X).*sin(pi*Y) ...
%                              + 13*D*pi^2*sin(2*pi*X).*sin(3*pi*Y) ...
%                              + 52*D*pi^2*sin(5*pi*X).*sin(pi*Y);
%
%                        3:  f =
%                     8*D*pi^2*sin(2*pi*X).*sin(3*pi*Y)*exp(-5*D*pi^2*t) ...
%                   - 24*D*pi^2*sin(5*pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t)
%
% - dx:                  grid spacing
% - dt:                  time step
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
% - X, Y:                grid points
% - timing_data:         array of the following form containing timing data
%                        [laplacian_construction time, solution time].
%                        If timing is not activated, timing_data is 
%                        set to [-1, -1].
% 
% NOTES:
% - The grid spacing is assumed to be the same in both the x and y directions.
%
% - The time step dt should be set to something less than dx/(16 D) so that
%   the temporal error does not dominate.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% CHANGE LOG:
% -----------
% 2007/09:  Initial version of code.
% 2008/02:  Changed code to be a function.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [u, u_exact, X, Y, timing_data] = ...
  solveDiffusionEqnCrankNicholson2d(D, ...
                                    source_term_type, ...
                                    dx, dt, ...
                                    t_final, ...
                                    debug_on, timing_on)

% check arguments
if (nargin < 5)
  error('solveDiffusionEqnCrankNicholson2d: missing arguments');
end
if (nargin < 6)
  debug_on = 0;
end
if (nargin < 7)
  timing_on = 0;
end

% construct grid
dy = dx;
N = 1/dx-1;
num_gridpts = N^2;
x = dx:dx:1-dx; x = x'; y = x;
[X,Y] = meshgrid(x,y);   % grid created with y as fastest direction
X = reshape(X,num_gridpts,1);
Y = reshape(Y,num_gridpts,1);

% set maximum number of grid cells to use for debug plotting
if (debug_on)
  N_plot = 101;
  if (N < N_plot)
    N_plot = N;
  end
  
  dx_plot = 2/(N_plot-1);
  x_plot = 0:dx_plot:1; y_plot = x_plot;
  [X_plot,Y_plot] = meshgrid(x_plot,y_plot);
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
Lx = kron(L_1D,speye(N));
Ly = kron(speye(N),L_1D);

% construct 5pt-Laplacian operator 
L = Lx + Ly; 

L_rhs = speye(N*N) + D*dt/2*L;
L_lhs = speye(N*N) - D*dt/2*L;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve diffusion equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% measure Laplacian construction time and
% restart clock to measure computation time
if (timing_on == 1)
  time_laplacian_construction = cputime - t_start;
  t_start = cputime;
end

% initialize concentration
disp('Initializing concentration field ...');
if (source_term_type == 1)
  u = 1 + Y/3 + X/2 + X.*Y/4 ...
    - sin(pi*X).*sin(pi*Y) + sin(2*pi*X).*sin(3*pi*Y);
elseif (source_term_type == 2)
  u = 1 + Y/3 + X/2 + X.*Y/4 + sin(5*pi*X).*sin(pi*Y);
elseif (source_term_type == 3) 
  u = 1 + Y/3 + X/2 + X.*Y/4 ...
    + sin(2*pi*X).*sin(3*pi*Y) - sin(5*pi*X).*sin(pi*Y);
else
  error('Invalid source term type.  Valid values are 1, 2, and 3.');
end


% Crank-Nicholson time integration
t = 0.0;
disp('Solving diffusion equation ...');
time_step = 1;
while (t < t_final)

  if (debug_on == 1)
  
    % compute exact solution
    if (source_term_type == 1)
      u_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
              - sin(pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t) ...
              + sin(2*pi*X).*sin(3*pi*Y)*exp(-13*D*pi^2*t);
    elseif (source_term_type == 2)
      u_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
              - sin(pi*X).*sin(pi*Y)*(1-exp(-2*D*pi^2*t)) ...
              + sin(2*pi*X).*sin(3*pi*Y)*(1-exp(-13*D*pi^2*t)) ...
              + sin(5*pi*X).*sin(pi*Y)*(2-exp(-26*D*pi^2*t));
    elseif (source_term_type == 3)
      u_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
              + sin(2*pi*X).*sin(3*pi*Y)*exp(-5*D*pi^2*t) ...
              - sin(5*pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t);
    else
      error('Invalid source term type');
    end
    err = u-u_exact;
    err_L_inf = norm(err,'inf');

    % add boundary conditions for plotting
    x_with_bdry = 0:dx:1; 
    x_with_bdry = x_with_bdry'; 
    y_with_bdry = x_with_bdry;   
    [X_with_bdry,Y_with_bdry] = meshgrid(x_with_bdry,y_with_bdry);
    N_with_bdry = N+2;
    u_with_bdry = zeros(N_with_bdry,N_with_bdry);
    u_with_bdry(1,1:N_with_bdry) = 1+x_with_bdry/2;  % y = 0
    u_with_bdry(N_with_bdry,1:N_with_bdry) = 4/3+3/4*x_with_bdry;  % y = 1
    u_with_bdry(1:N_with_bdry,1) = 1+y_with_bdry/3;  % x = 0
    u_with_bdry(1:N_with_bdry,N_with_bdry) = 1.5+7/12*y_with_bdry;  % x = 1
    u_with_bdry(2:N_with_bdry-1,2:N_with_bdry-1) = reshape(u,N,N);

    err_with_bdry = zeros(N_with_bdry,N_with_bdry);
    err_with_bdry(2:N_with_bdry-1,2:N_with_bdry-1) = reshape(err,N,N);


    % plot current solution
    figure(1); clf;
    u_plot = interp2(reshape(X_with_bdry,N_with_bdry,N_with_bdry), ...
                     reshape(Y_with_bdry,N_with_bdry,N_with_bdry), ...
                     reshape(u_with_bdry,N_with_bdry,N_with_bdry), ...
                     X_plot,Y_plot,'*cubic');
    surf(x_plot,y_plot,u_plot);
    title_string = sprintf('t = %f',t);
    title(title_string);
    xlabel('x'); ylabel('y');

    % plot current error
    figure(2); clf;
    err_plot = interp2(reshape(X_with_bdry,N_with_bdry,N_with_bdry), ...
                       reshape(Y_with_bdry,N_with_bdry,N_with_bdry), ...
                       reshape(err_with_bdry,N_with_bdry,N_with_bdry), ...
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
    L_rhs = speye(N*N) + D*dt/2*L;
    L_lhs = speye(N*N) - D*dt/2*L;
  end

  % compute source terms
  if (source_term_type == 1)
    f = 0;
  elseif (source_term_type == 2)
    f = -2*D*pi^2*sin(pi*X).*sin(pi*Y) ...
      + 13*D*pi^2*sin(2*pi*X).*sin(3*pi*Y) ...
      + 52*D*pi^2*sin(5*pi*X).*sin(pi*Y);
  elseif (source_term_type == 3)
    f_cur = 8*D*pi^2*sin(2*pi*X).*sin(3*pi*Y)*exp(-5*D*pi^2*t) ...
          - 24*D*pi^2*sin(5*pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t);
    f_next = 8*D*pi^2*sin(2*pi*X).*sin(3*pi*Y)*exp(-5*D*pi^2*(t+dt)) ...
           - 24*D*pi^2*sin(5*pi*X).*sin(pi*Y)*exp(-2*D*pi^2*(t+dt));
    f = 0.5*(f_cur+f_next);
  else
    error('Invalid source term type');
  end

  % compute RHS
  rhs = L_rhs*u + dt*f;

  % impose boundary conditions
  rhs(1:N:num_gridpts) = rhs(1:N:num_gridpts) + D*dt/dy/dy*(1+x/2);  % y = 0
  rhs(N:N:num_gridpts) = ...
    rhs(N:N:num_gridpts) + D*dt/dy/dy*(4/3+3/4*x);  % y = 1
  rhs(1:N) = rhs(1:N) + D*dt/dx/dx*(1+y/3);  % x = 0
  rhs(num_gridpts-N+1:num_gridpts) = ...
    rhs(num_gridpts-N+1:num_gridpts) + D*dt/dx/dx*(1.5+7/12*y);  % x = 1

  % update solution
  u = L_lhs\rhs;

  % update time and time step
  t = t + dt;
  time_step = time_step + 1;

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

% add boundaries to grid
x_with_bdry = 0:dx:1; 
x_with_bdry = x_with_bdry'; 
y_with_bdry = x_with_bdry;   
[X,Y] = meshgrid(x_with_bdry,y_with_bdry);
N_with_bdry = N+2;
X = reshape(X, N_with_bdry^2, 1);
Y = reshape(Y, N_with_bdry^2, 1);

% compute exact solution
if (source_term_type == 1)
  u_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
          - sin(pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t) ...
          + sin(2*pi*X).*sin(3*pi*Y)*exp(-13*D*pi^2*t);
elseif (source_term_type == 2)
  u_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
          - sin(pi*X).*sin(pi*Y)*(1-exp(-2*D*pi^2*t)) ...
          + sin(2*pi*X).*sin(3*pi*Y)*(1-exp(-13*D*pi^2*t)) ...
          + sin(5*pi*X).*sin(pi*Y)*(2-exp(-26*D*pi^2*t));
elseif (source_term_type == 3)
  u_exact = 1 + Y/3 + X/2 + X.*Y/4 ...
          + sin(2*pi*X).*sin(3*pi*Y)*exp(-5*D*pi^2*t) ...
          - sin(5*pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t);
else
  error('Invalid source term type');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add boundary conditions on numerical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_with_bdry = zeros(N_with_bdry,N_with_bdry);
u_with_bdry(1,1:N_with_bdry) = 1+x_with_bdry/2;  % y = 0
u_with_bdry(N_with_bdry,1:N_with_bdry) = 4/3+3/4*x_with_bdry;  % y = 1
u_with_bdry(1:N_with_bdry,1) = 1+y_with_bdry/3;  % x = 0
u_with_bdry(1:N_with_bdry,N_with_bdry) = 1.5+7/12*y_with_bdry;  % x = 1
u_with_bdry(2:N_with_bdry-1,2:N_with_bdry-1) = reshape(u,N,N);
u = reshape(u_with_bdry, N_with_bdry^2, 1);

