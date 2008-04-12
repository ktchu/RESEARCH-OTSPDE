%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqnForwardEulerOTS2d() computes the solutions of the 2d
% diffusion equation
%
%   u_t = D (u_xx + u_yy) + f
%
% on an irregular domain bounded by the zero level set of a user-specified
% function phi with Dirichlet boundary conditions imposed everywhere on the 
% boundary.  The computational domain is taken to the square -1 < x,y, 1. 
% The numerical solution is computed on a node-centered grid using 
% forward Euler time integration with the standard second-order 9pt 
% discretization of the Laplacian.  A higher-order correction for the
% source term is used.  The optimal time step is dt = dx^2/(6D) and
% correction terms are used.
%
% The Dirichlet boundary conditions are provided by the exact analytical
% solution to the problem on the domain -1 < x < 1.  The value of u
% in the ghostcells may be set to linear, quadratic, or cubic via the
% variable bc_interpolant_order.
%
% USAGE:
%  function [u, u_exact, X, Y, timing_data, idx_ghostcells], ...
%    solveDiffusionEqnForwardEulerOTS2d(D, ...
%                                       source_term_type, ...
%                                       phi, ...
%                                       dx, ...
%                                       t_final, ...
%                                       bc_interpolant_order, ...
%                                       debug_on, timing_on, ...
%                                       zero_level_set_tol, extrap_tol)
%
% Arguments:
% - D:                     diffusion coefficient
% - source_term_type:      type of source term to use in computation.
%                          1:  f = 0
%
%                          2:  f = -2*D*pi^2*sin(pi*X).*sin(pi*Y) ...
%                                + 13*D*pi^2*sin(2*pi*X).*sin(3*pi*Y) ...
%                                + 52*D*pi^2*sin(5*pi*X).*sin(pi*Y);
%
%                          3:  f =
%                     8*D*pi^2*sin(2*pi*X).*sin(3*pi*Y)*exp(-5*D*pi^2*t) ...
%                   - 24*D*pi^2*sin(5*pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t)
%
% - phi:                   function whose zero level set defines the 
%                          boundary of the domain.  phi may be in either
%                          grid or vector form.
% - dx:                    grid spacing
% - t_final:               final time
% - bc_interpolant_order:  order of polynomial interpolant used to 
%                          extend boundary conditions from zero level 
%                          set to ghost cells.  Linear (=1), quadratic (=2) 
%                          and cubic (=3) interpolation are supported.
%                          (default = 3)
% - debug_on:              flag indicating whether debugging information
%                          should be displayed.  To turn on debugging,
%                          set debug_on to 1.
%                          (default = 0)
% - timing_on:             flag indicating whether timing information
%                          should be collected.  To activate timing,
%                          set timing_on to 1.
%                          (default = 0)
% - zero_level_set_tol:    tolerance used to identify grid points that
%                          should be treated as being on the zero level set.
%                          If no tolerance is specified, the default value
%                          set by computeIrregBdryParams() is used.  To use 
%                          the default value, set zero_level_set_tol to -1.
% - extrap_tol:            tolerance used to deside when the grid points
%                          used for extrapolation should be shifted towards
%                          the interior of the domain.  The shift ensures
%                          that the Lagrange interpolation formula used to
%                          fill edge ghost cells is stable.  If no tolerance 
%                          is specified, the default value set by 
%                          computeIrregBdryParams() is used.  To use the
%                          default value, set extrap_tol to -1.
%
% Return values:
% - u:                   numerical solution
% - u_exact:             analytical solution
% - X, Y:                x- and y-coordinates of grid points
% - timing_data:         array of the following form containing timing data
%                        [BC parameter calculation time, laplacian
%                         construction time, solution time].
%                        set to [-1, -1, -1].
% - idx_ghostcells:      cell array containing the indices of the 
%                        grid points on the boundary, edge ghost cells,
%                        and corner ghostcells
%
% NOTES:
% - The grid spacing is assumed to be the same in both the x and y directions.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



function [u, u_exact, X, Y, timing_data, idx_ghostcells] = ...
  solveDiffusionEqnForwardEulerOTS2d(D, ...
                                     source_term_type, ...
                                     phi, ...
                                     dx, ...
                                     t_final, ...
                                     bc_interpolant_order, ...
                                     debug_on, timing_on, ...
                                     zero_level_set_tol, extrap_tol)

% check arguments
if (nargin < 5)
  error('solveDiffusionEqnForwardEulerOTS2d: missing arguments');
end
if (nargin < 6)
  bc_interpolant_order = 3;
end
if (nargin < 7)
  debug_on = 0;
end
if (nargin < 8)
  timing_on = 0;
end
if (nargin < 9)
  zero_level_set_tol = -1;
end
if (nargin < 10)
  extrap_tol = -1;
end

% construct grid
dy = dx;
N = 2.0/dx + 1;
num_gridpts = N*N;
x = -1:dx:1; x = x'; y = x;
[X,Y] = meshgrid(x,y);   % grid created with y as fastest direction
X = reshape(X,num_gridpts,1);
Y = reshape(Y,num_gridpts,1);

% compute optimal time step 
dt = dx^2/D/6;

% convert phi to a vector if it is in matrix form
if (size(phi,1) ~=1 & size(phi,2) ~=1)
  phi = reshape(phi,num_gridpts,1);
end

% set maximum number of grid cells to use for debug plotting
if (debug_on)
  N_plot = 101;
  if (N < N_plot)
    N_plot = N;
  end

  dx_plot = 2/(N_plot-1);
  x_plot = -1:dx_plot:1; y_plot = x_plot;
  [X_plot,Y_plot] = meshgrid(x_plot,y_plot);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for main computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start clock to measure time to compute boundary condition parameters
disp('Calculating parameters for boundary conditions at irregular grid points ...');
if (timing_on == 1)
  t_start = cputime;
end

% compute parameters for irregular boundary grid points
[idx_inside, idx_outside, idx_bdry, ...
 idx_ghostcells, idx_ghostcells_edge, idx_ghostcells_corner, ...
 idx_edge_bc_interp_pts, coefs_edge_bc_interp_pts, coefs_edge_bc_bdry_pts, ...
 x_edge_bc_bdry_pts, y_edge_bc_bdry_pts, ...
 idx_corner_bc_interp_pts, coefs_corner_bc_interp_pts] = ...    
  computeIrregBdryParams(X, Y, N, dx, phi, ...
                         bc_interpolant_order, ...
                         zero_level_set_tol, extrap_tol);

% measure time to calculate boundary condition parameters and 
% restart clock to measure Laplacian construction time 
if (timing_on == 1)
  time_bc_params = cputime - t_start;
  t_start = cputime;
end

% construct Laplacian operators (with boundary conditions) 
disp('Constructing Laplacian operators ...');
e = ones(N,1);
L_1D = 1/dx^2*spdiags([e -2*e e], -1:1, N, N);
Lx = kron(L_1D,speye(N));
Ly = kron(speye(N),L_1D);

% construct 5pt-Laplacian operator (with boundary conditions) 
L_5pt = Lx + Ly; 
L_5pt(1:N:N*N,:) = 0;         % Dirichlet BC x = 0
L_5pt(N:N:N*N,:) = 0;         % Dirichlet BC x = 1
L_5pt(1:N,:) = 0;             % Dirichlet BC y = 0
L_5pt((N-1)*N+1:N*N,:) = 0;   % Dirichlet BC y = 1

% construct 9pt-Laplacian operator (with boundary conditions) 
L_9pt = L_5pt + dx*dx/6*Lx*Ly;
L_9pt(1:N:N*N,:) = 0;         % Dirichlet BC x = 0
L_9pt(N:N:N*N,:) = 0;         % Dirichlet BC x = 1
L_9pt(1:N,:) = 0;             % Dirichlet BC y = 0
L_9pt((N-1)*N+1:N*N,:) = 0;   % Dirichlet BC y = 1


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
  error('Invalid source term type');
end
u(idx_outside) = 0;

% forward Euler time integration
disp('Solving diffusion equation ...');
t = 0.0;
time_step = 1;
while (t < t_final)
 
  if (debug_on == 1)
 
    % impose boundary conditions
    x_on_bdry = X(idx_bdry); y_on_bdry = Y(idx_bdry);
    if (source_term_type == 1)
      u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
             + x_on_bdry.*y_on_bdry/4 ...
             - sin(pi*x_on_bdry).*sin(pi*y_on_bdry) ...
              *exp(-2*D*pi^2*t) ...
             + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
              *exp(-13*D*pi^2*t); 
    elseif (source_term_type == 2)
      u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
             + x_on_bdry.*y_on_bdry/4 ...
             - sin(pi*x_on_bdry).*sin(pi*y_on_bdry) ...
              *(1-exp(-2*D*pi^2*t)) ...
             + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
              *(1-exp(-13*D*pi^2*t)) ...
             + sin(5*pi*x_on_bdry).*sin(pi*y_on_bdry) ...
              *(2-exp(-26*D*pi^2*t));
    elseif (source_term_type == 3)
      u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
             + x_on_bdry.*y_on_bdry/4 ...
             + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
              *exp(-5*D*pi^2*t) ...
             - sin(5*pi*x_on_bdry).*sin(pi*y_on_bdry) ...
              *exp(-2*D*pi^2*t);
    else
      error('Invalid source term type');
    end
    u(idx_bdry) = u_bdry;

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
    err = abs(u-u_exact);
    err(idx_outside) = 0;
    err_L_inf = norm(err,'inf');
    err_L2    = norm(err,2)*dx;

    % plot current solution
    figure(1); clf;
    u(idx_outside) = 0;
    u_plot = interp2(reshape(X,N,N),reshape(Y,N,N),reshape(u,N,N), ...
                     X_plot,Y_plot,'*cubic');
    surf(x_plot,y_plot,u_plot);
    title_string = sprintf('t = %f',t);
    title(title_string);
    xlabel('x'); ylabel('y');

    % plot current error
    figure(2); clf;
    err_plot = interp2(reshape(X,N,N),reshape(Y,N,N),reshape(err,N,N), ...
                       X_plot,Y_plot,'*cubic');
    surf(x_plot,y_plot,err_plot);
    title_string = sprintf('t = %f',t);
    title(title_string);
    xlabel('x'); ylabel('y');
    drawnow

    % display current status of solution
    disp('------------------------------------------------------------------');
    status_str1 = sprintf('Time = %f, Time step = %d', t, time_step);
    status_str2 = sprintf('L inf Error in Solution = %g', err_L_inf);
    status_str3 = sprintf('L2 Error in Solution = %g', err_L2);
    disp(status_str1);
    disp(status_str2);
    disp(status_str3);

  end

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t;
  end

  % impose boundary conditions
  x_on_bdry = X(idx_bdry); y_on_bdry = Y(idx_bdry);
  if (source_term_type == 1)
    u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
           + x_on_bdry.*y_on_bdry/4 ...
           - sin(pi*x_on_bdry).*sin(pi*y_on_bdry) ...
            *exp(-2*D*pi^2*t) ...
           + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
            *exp(-13*D*pi^2*t); 
  elseif (source_term_type == 2)
    u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
           + x_on_bdry.*y_on_bdry/4 ...
           - sin(pi*x_on_bdry).*sin(pi*y_on_bdry) ...
            *(1-exp(-2*D*pi^2*t)) ...
           + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
            *(1-exp(-13*D*pi^2*t)) ...
           + sin(5*pi*x_on_bdry).*sin(pi*y_on_bdry) ...
            *(2-exp(-26*D*pi^2*t));
  elseif (source_term_type == 3)
    u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
           + x_on_bdry.*y_on_bdry/4 ...
           + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
            *exp(-5*D*pi^2*t) ...
           - sin(5*pi*x_on_bdry).*sin(pi*y_on_bdry) ...
            *exp(-2*D*pi^2*t);
  else
    error('Invalid source term type');
  end

  if (source_term_type == 1)
    u_bdry_edge = 1 + y_edge_bc_bdry_pts/3 + x_edge_bc_bdry_pts/2 ...
                + x_edge_bc_bdry_pts.*y_edge_bc_bdry_pts/4 ...
                - sin(pi*x_edge_bc_bdry_pts).*sin(pi*y_edge_bc_bdry_pts) ...
                 *exp(-2*D*pi^2*t) ...
                + sin(2*pi*x_edge_bc_bdry_pts).*sin(3*pi*y_edge_bc_bdry_pts) ...
                 *exp(-13*D*pi^2*t); 
  elseif (source_term_type == 2)
    u_bdry_edge = 1 + y_edge_bc_bdry_pts/3 + x_edge_bc_bdry_pts/2 ...
                + x_edge_bc_bdry_pts.*y_edge_bc_bdry_pts/4 ...
                - sin(pi*x_edge_bc_bdry_pts).*sin(pi*y_edge_bc_bdry_pts) ...
                 *(1-exp(-2*D*pi^2*t)) ...
                + sin(2*pi*x_edge_bc_bdry_pts).*sin(3*pi*y_edge_bc_bdry_pts) ...
                 *(1-exp(-13*D*pi^2*t)) ...
                + sin(5*pi*x_edge_bc_bdry_pts).*sin(pi*y_edge_bc_bdry_pts) ...
                 *(2-exp(-26*D*pi^2*t));
  elseif (source_term_type == 3)
    u_bdry_edge = 1 + y_edge_bc_bdry_pts/3 + x_edge_bc_bdry_pts/2 ...
                + x_edge_bc_bdry_pts.*y_edge_bc_bdry_pts/4 ...
                + sin(2*pi*x_edge_bc_bdry_pts).*sin(3*pi*y_edge_bc_bdry_pts) ...
                 *exp(-5*D*pi^2*t) ...
                - sin(5*pi*x_edge_bc_bdry_pts).*sin(pi*y_edge_bc_bdry_pts) ...
                 *exp(-2*D*pi^2*t);
  else
    error('Invalid source term type');
  end

  % set concentration of points on the boundary
  u(idx_bdry) = u_bdry;

  % fill edge ghost cells by extrapolating from interior
  if ( bc_interpolant_order == 1 )

    u(idx_ghostcells_edge) = ...
        coefs_edge_bc_bdry_pts.*u_bdry_edge ...
      + coefs_edge_bc_interp_pts.*u(idx_edge_bc_interp_pts);

  elseif ( bc_interpolant_order == 2 )

    u(idx_ghostcells_edge) = ...
        coefs_edge_bc_bdry_pts.*u_bdry_edge ...
      + coefs_edge_bc_interp_pts(:,1).*u(idx_edge_bc_interp_pts(:,1)) ...
      + coefs_edge_bc_interp_pts(:,2).*u(idx_edge_bc_interp_pts(:,2));

  elseif ( bc_interpolant_order == 3 )

    u(idx_ghostcells_edge) = ...
        coefs_edge_bc_bdry_pts.*u_bdry_edge ...
      + coefs_edge_bc_interp_pts(:,1).*u(idx_edge_bc_interp_pts(:,1)) ...
      + coefs_edge_bc_interp_pts(:,2).*u(idx_edge_bc_interp_pts(:,2)) ...
      + coefs_edge_bc_interp_pts(:,3).*u(idx_edge_bc_interp_pts(:,3));

  else

    error('Invalid order for polynomial interpolant.');

  end

  % fill corner ghost cells by extrapolating from interior
  % and edge ghost cells
  if ( bc_interpolant_order == 1 ...
     | bc_interpolant_order == 2 )

    % linear extrapolation for corner points is sufficient
    % to obtain a second- or third-order accurate solution

    u(idx_ghostcells_corner) = ...
        coefs_corner_bc_interp_pts(:,1).*u(idx_corner_bc_interp_pts(:,1)) ...
      + coefs_corner_bc_interp_pts(:,2).*u(idx_corner_bc_interp_pts(:,2)) ...
      + coefs_corner_bc_interp_pts(:,3).*u(idx_corner_bc_interp_pts(:,3));

  elseif ( bc_interpolant_order == 3 )

    % quadratic extrapolation for corner points is required
    % to obtain a fourth-order accurate solution

    u(idx_ghostcells_corner) = ...
        coefs_corner_bc_interp_pts(:,1).*u(idx_corner_bc_interp_pts(:,1)) ...
      + coefs_corner_bc_interp_pts(:,2).*u(idx_corner_bc_interp_pts(:,2)) ...
      + coefs_corner_bc_interp_pts(:,3).*u(idx_corner_bc_interp_pts(:,3)) ...
      + coefs_corner_bc_interp_pts(:,4).*u(idx_corner_bc_interp_pts(:,4)) ...
      + coefs_corner_bc_interp_pts(:,5).*u(idx_corner_bc_interp_pts(:,5)) ...
      + coefs_corner_bc_interp_pts(:,6).*u(idx_corner_bc_interp_pts(:,6)) ...
      + coefs_corner_bc_interp_pts(:,7).*u(idx_corner_bc_interp_pts(:,7)) ...
      + coefs_corner_bc_interp_pts(:,8).*u(idx_corner_bc_interp_pts(:,8));

  else

    error('Invalid order for polynomial interpolant.');

  end

  % compute error in ghostcells 
  if (debug_on == 1)
    err_bdry_ghostcells = ...
      u(idx_bdry)-u_exact(idx_bdry);
    err_L_inf_bdry_ghostcells = norm(err_bdry_ghostcells,inf);
    err_edge_ghostcells = ...
      u(idx_ghostcells_edge)-u_exact(idx_ghostcells_edge);
    err_L_inf_edge_ghostcells = norm(err_edge_ghostcells,inf);
    err_corner_ghostcells = ...
      u(idx_ghostcells_corner)-u_exact(idx_ghostcells_corner);
    err_L_inf_corner_ghostcells = norm(err_corner_ghostcells,inf);
    bdry_bc_err_str = ...
      sprintf('L inf Error in Bdry Ghostcells = %g', ...
              err_L_inf_bdry_ghostcells);
    edge_bc_err_str = ...
      sprintf('L inf Error in Edge Ghostcells = %g', ...
              err_L_inf_edge_ghostcells);
    corner_bc_err_str = ...
      sprintf('L inf Error in Corner Ghostcells = %g', ...
              err_L_inf_corner_ghostcells);
    disp(bdry_bc_err_str);
    disp(edge_bc_err_str);
    disp(corner_bc_err_str);
  end

  % compute source terms
  if (source_term_type == 1)
    f = 0;
    f_t = 0;
  elseif (source_term_type == 2)
    f = -2*D*pi^2*sin(pi*X).*sin(pi*Y) ...
      + 13*D*pi^2*sin(2*pi*X).*sin(3*pi*Y) ...
      + 52*D*pi^2*sin(5*pi*X).*sin(pi*Y);
    f_t = 0;
  elseif (source_term_type == 3)
    f = 8*D*pi^2*sin(2*pi*X).*sin(3*pi*Y)*exp(-5*D*pi^2*t) ...
      - 24*D*pi^2*sin(5*pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t);
    f_t = -40*D^2*pi^4*sin(2*pi*X).*sin(3*pi*Y)*exp(-5*D*pi^2*t) ...
        + 48*D^2*pi^4*sin(5*pi*X).*sin(pi*Y)*exp(-2*D*pi^2*t);
  else
    error('Invalid source term type');
  end

  % update solution
  if (source_term_type == 1)
    u = u + dt*D*(L_9pt*u);
  elseif (source_term_type == 2)
    u = u + dt*D*(L_9pt*u) + dt*f + 0.5*dt^2*D*L_5pt*f; 
  elseif (source_term_type == 3)
    u = u + dt*D*(L_9pt*u) + dt*f + 0.5*dt^2*(f_t + D*L_5pt*f); 
  else 
    error('Invalid source term selection');
  end
  u(idx_outside) = 0;

  % update time and time_step
  t = t + dt;
  time_step = time_step + 1;

end

% measure time to solve diffusion equation
if (timing_on == 1)
  time_solve = cputime - t_start;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output timing statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (timing_on)
  timing_data = [time_bc_params, time_laplacian_construction, time_solve];

  if (debug_on)
    disp(' ');
    disp('==================================================================');
    disp('Computation Statistics');
    if (timing_on == 1)
      bc_param_time_str = ...
        sprintf('  Boundary Condition Parameter Calculation Time: %f', ...
        time_bc_params);
      lapl_construct_time_str = ...
        sprintf('  Laplacian Construction Time: %f', time_laplacian_construction);
      comp_time_str = sprintf('  Solution Time: %f', time_solve);
      disp(bc_param_time_str);
      disp(lapl_construct_time_str);
      disp(comp_time_str);
    end
    disp('==================================================================');

  end

else
  
  timing_data = [-1, -1, -1];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set boundary conditions on numerical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% impose boundary conditions
x_on_bdry = X(idx_bdry); y_on_bdry = Y(idx_bdry);
if (source_term_type == 1)
  u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
         + x_on_bdry.*y_on_bdry/4 ...
         - sin(pi*x_on_bdry).*sin(pi*y_on_bdry) ...
          *exp(-2*D*pi^2*t) ...
         + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
          *exp(-13*D*pi^2*t); 
elseif (source_term_type == 2)
  u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
         + x_on_bdry.*y_on_bdry/4 ...
         - sin(pi*x_on_bdry).*sin(pi*y_on_bdry) ...
          *(1-exp(-2*D*pi^2*t)) ...
         + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
          *(1-exp(-13*D*pi^2*t)) ...
         + sin(5*pi*x_on_bdry).*sin(pi*y_on_bdry) ...
          *(2-exp(-26*D*pi^2*t));
elseif (source_term_type == 3)
  u_bdry = 1 + y_on_bdry/3 + x_on_bdry/2 ...
         + x_on_bdry.*y_on_bdry/4 ...
         + sin(2*pi*x_on_bdry).*sin(3*pi*y_on_bdry) ...
          *exp(-5*D*pi^2*t) ...
         - sin(5*pi*x_on_bdry).*sin(pi*y_on_bdry) ...
          *exp(-2*D*pi^2*t);
else
  error('Invalid source term type');
end
u(idx_bdry) = u_bdry;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute exact solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
u_exact(idx_outside) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set idx_ghostcells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_ghostcells = {idx_bdry, idx_ghostcells_edge, idx_ghostcells_corner};
