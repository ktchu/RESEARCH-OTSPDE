%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script computes the solutions of the 1d heat equation 
%
%   T_t = kappa T_xx + f
%
% on the domain x_lo < x < x_hi with a Dirichlet boundary condition imposed 
% at both boundaries where 0 < x_lo, x_hi < 1.  The computational domain is 
% taken to be 0 < x < 1 and the mathematical boundaries of the domain are 
% treated as irregular grid points.  The numerical solution is computed on 
% a node-centered grid using forward Euler time integration with a 
% second-order central difference approximation to the Laplacian.  A 
% higher-order correction for the source term is used.  The optimal time 
% step is dt = dx^2/kappa/6.
%  
% The Dirichlet boundary conditions are provided by the exact analytical
% solution to the problem on the domain 0 < x < 1.  The value of T
% in the ghostcells may be set to constant, linear, quadratic, or cubic 
% via the variable bc_polynomial_interpolant_order.
%
% There are user-specified two choices for the source term f:
%
% (1) f = 0
% (2) f = kappa*(50*sin(3/2*pi*x) - 70*sin(15/2*pi*x) + 100*sin(21/2*pi*x))
%  
% Kevin T. Chu
% 2007 September
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% Set timing_run to 1 to disable error calculation and plotting
% during time integration loop.
timing_run = 0;

% source term
source_term_type = 2;

% order of polynomial interpolant for imposing boundary condition 
bc_polynomial_interpolant_order = 3;

% physical parameters
kappa = 2.0;  % thermal conductivity

% time integration parameters
t_init  = 0.0;
t_final = 0.001;

% number of grid cells to use to discretize the full computational 0 < x < 1
% NOTE: number of grid points is N = (N_cell+1)
N_cell = 50;
N_cell = 200;
N_cell = 100;

% boundary points
x_lo = 1/12;
x_hi = 5/6;

% construct grid
N = N_cell+1;
dx = 1/N_cell;
x = 0:dx:1;  x = x';   

% set dt
dt = dx^2/6/kappa;

% compute indices of 
% (1) interior and exterior grid points
% (2) ghostcells
x_bdry = [x_lo; x_hi];
idx_inside  = find(x > x_lo & x < x_hi);
idx_outside = find(x <= x_lo | x >= x_hi);
idx_ghostcells = find( (abs(x-x_lo) < dx & x < x_lo) ...
                     | (abs(x-x_hi) < dx & x > x_hi) );

% precompute several quantities required for extrapolation of boundary
% data to ghostcells:
% (1) indices of interior neighbor points used for extrapolation,
% (2) coefficients for extrapolation
% (3) coefficients for extrapolation
extrap_tol = dx^2;
if (bc_polynomial_interpolant_order == 0)
  idx_bc_interp_pts   = [];
  coefs_bc_interp_pts = [];
  coefs_bc_bdry_pts   = ones(2,1);
elseif (bc_polynomial_interpolant_order == 1)

  idx_bc_interp_pts   = zeros(2,1);
  coefs_bc_interp_pts = zeros(2,1);
  coefs_bc_bdry_pts   = zeros(2,1);

  % x = x_lo
  x_b = x_bdry(1);
  x_g = x(idx_ghostcells(1));
  idx_bc_interp_pts(1) = idx_ghostcells(1)+1;
  x_1 = x(idx_bc_interp_pts(1));
  if ( x_1 - x_b < extrap_tol)
    idx_bc_interp_pts(1) = idx_ghostcells(1)+2;
    x_1 = x(idx_bc_interp_pts(1));
  end
  coefs_bc_bdry_pts(1)   = (x_g-x_1)/(x_b-x_1);
  coefs_bc_interp_pts(1) = (x_g-x_b)/(x_1-x_b);

  % x = x_hi
  x_b = x_bdry(2);
  x_g = x(idx_ghostcells(2));
  idx_bc_interp_pts(2) = idx_ghostcells(2)-1;
  x_1 = x(idx_bc_interp_pts(2));
  if (x_b - x_1 < extrap_tol)
    idx_bc_interp_pts(2) = idx_ghostcells(2)-2;
    x_1 = x(idx_bc_interp_pts(2));
  end
  coefs_bc_bdry_pts(2)   = (x_g-x_1)/(x_b-x_1);
  coefs_bc_interp_pts(2) = (x_g-x_b)/(x_1-x_b);

elseif (bc_polynomial_interpolant_order == 2)

  idx_bc_interp_pts   = zeros(2,2);
  coefs_bc_interp_pts = zeros(2,2);
  coefs_bc_bdry_pts   = zeros(2,1);

  % x = x_lo
  x_b = x_bdry(1);
  x_g = x(idx_ghostcells(1));
  idx_bc_interp_pts(1,1) = idx_ghostcells(1)+1;
  idx_bc_interp_pts(1,2) = idx_ghostcells(1)+2;
  x_1 = x(idx_bc_interp_pts(1,1));
  x_2 = x(idx_bc_interp_pts(1,2));
  if (x_1 - x_b < extrap_tol)
    idx_bc_interp_pts(1,1) = idx_ghostcells(1)+2;
    idx_bc_interp_pts(1,2) = idx_ghostcells(1)+3;
    x_1 = x(idx_bc_interp_pts(1,1));
    x_2 = x(idx_bc_interp_pts(1,2));
  end
  coefs_bc_bdry_pts(1)     = (x_g-x_1)*(x_g-x_2)/(x_b-x_1)/(x_b-x_2);
  coefs_bc_interp_pts(1,1) = (x_g-x_b)*(x_g-x_2)/(x_1-x_b)/(x_1-x_2);
  coefs_bc_interp_pts(1,2) = (x_g-x_b)*(x_g-x_1)/(x_2-x_b)/(x_2-x_1);

  % x = x_hi
  x_b = x_bdry(2);
  x_g = x(idx_ghostcells(2));
  idx_bc_interp_pts(2,1) = idx_ghostcells(2)-1;
  idx_bc_interp_pts(2,2) = idx_ghostcells(2)-2;
  x_1 = x(idx_bc_interp_pts(2,1));
  x_2 = x(idx_bc_interp_pts(2,1));
  if (x_b - x_1 < extrap_tol)
    idx_bc_interp_pts(2,1) = idx_ghostcells(2)-2;
    idx_bc_interp_pts(2,2) = idx_ghostcells(2)-3;
    x_1 = x(idx_bc_interp_pts(2,1));
    x_2 = x(idx_bc_interp_pts(2,2));
  end
  coefs_bc_bdry_pts(2)     = (x_g-x_1)*(x_g-x_2)/(x_b-x_1)/(x_b-x_2);
  coefs_bc_interp_pts(2,1) = (x_g-x_b)*(x_g-x_2)/(x_1-x_b)/(x_1-x_2);
  coefs_bc_interp_pts(2,2) = (x_g-x_b)*(x_g-x_1)/(x_2-x_b)/(x_2-x_1);

elseif (bc_polynomial_interpolant_order == 3)

  idx_bc_interp_pts   = zeros(2,3);
  coefs_bc_interp_pts = zeros(2,3);
  coefs_bc_bdry_pts   = zeros(2,1);

  % x = x_lo
  x_b = x_bdry(1);
  x_g = x(idx_ghostcells(1));
  idx_bc_interp_pts(1,1) = idx_ghostcells(1)+1;
  idx_bc_interp_pts(1,2) = idx_ghostcells(1)+2;
  idx_bc_interp_pts(1,3) = idx_ghostcells(1)+3;
  x_1 = x(idx_bc_interp_pts(1,1));
  x_2 = x(idx_bc_interp_pts(1,2));
  x_3 = x(idx_bc_interp_pts(1,3));
  if (x_1 - x_b < extrap_tol)
    idx_bc_interp_pts(1,1) = idx_ghostcells(1)+2;
    idx_bc_interp_pts(1,2) = idx_ghostcells(1)+3;
    idx_bc_interp_pts(1,3) = idx_ghostcells(1)+4
    x_1 = x(idx_bc_interp_pts(1,1));
    x_2 = x(idx_bc_interp_pts(1,2));
    x_3 = x(idx_bc_interp_pts(1,3));
  end
  coefs_bc_bdry_pts(1)     = (x_g-x_1)*(x_g-x_2)*(x_g-x_3)...
                           / (x_b-x_1)/(x_b-x_2)/(x_b-x_3);
  coefs_bc_interp_pts(1,1) = (x_g-x_b)*(x_g-x_2)*(x_g-x_3)...
                           / (x_1-x_b)/(x_1-x_2)/(x_1-x_3);
  coefs_bc_interp_pts(1,2) = (x_g-x_b)*(x_g-x_1)*(x_g-x_3)...
                           / (x_2-x_b)/(x_2-x_1)/(x_2-x_3);
  coefs_bc_interp_pts(1,3) = (x_g-x_b)*(x_g-x_1)*(x_g-x_2)...
                           / (x_3-x_b)/(x_3-x_1)/(x_3-x_2);

  % x = x_hi
  x_b = x_bdry(2);
  x_g = x(idx_ghostcells(2));
  idx_bc_interp_pts(2,1) = idx_ghostcells(2)-1;
  idx_bc_interp_pts(2,2) = idx_ghostcells(2)-2;
  idx_bc_interp_pts(2,3) = idx_ghostcells(2)-3;
  x_1 = x(idx_bc_interp_pts(2,1));
  x_2 = x(idx_bc_interp_pts(2,2));
  x_3 = x(idx_bc_interp_pts(2,3));
  if (x_b - x_1 < extrap_tol)
    idx_bc_interp_pts(2,1) = idx_ghostcells(2)-2;
    idx_bc_interp_pts(2,2) = idx_ghostcells(2)-3;
    idx_bc_interp_pts(2,3) = idx_ghostcells(2)-4;
    x_1 = x(idx_bc_interp_pts(2,1));
    x_2 = x(idx_bc_interp_pts(2,2));
    x_3 = x(idx_bc_interp_pts(2,3));
  end
  coefs_bc_bdry_pts(2)     = (x_g-x_1)*(x_g-x_2)*(x_g-x_3)...
                           / (x_b-x_1)/(x_b-x_2)/(x_b-x_3);
  coefs_bc_interp_pts(2,1) = (x_g-x_b)*(x_g-x_2)*(x_g-x_3)...
                           / (x_1-x_b)/(x_1-x_2)/(x_1-x_3);
  coefs_bc_interp_pts(2,2) = (x_g-x_b)*(x_g-x_1)*(x_g-x_3)...
                           / (x_2-x_b)/(x_2-x_1)/(x_2-x_3);
  coefs_bc_interp_pts(2,3) = (x_g-x_b)*(x_g-x_1)*(x_g-x_2)...
                           / (x_3-x_b)/(x_3-x_1)/(x_3-x_2);

else
  error('Invalid order for polynomial interpolant.');
end

% construct Laplacian operator (with boundary conditions) 
L = (diag(ones(N-1,1),1) - 2*diag(ones(N,1),0) + diag(ones(N-1,1),-1))/dx/dx;
L(idx_outside,:) = 0; % discard part of Laplacian for exterior points

% initialize temperature
T = 1.0 + 0.5*x + 2*sin(5/2*pi*x) - 4*sin(11/2*pi*x) + 3*sin(7/2*pi*x);
T(idx_outside) = 0;

% compute source terms
if (source_term_type == 1)
  f = 0;
elseif (source_term_type == 2)
  f = kappa*(50*sin(3/2*pi*x) - 70*sin(15/2*pi*x) + 100*sin(21/2*pi*x));
else
  error('Invalid source term type.');
end

% forward Euler time integration
t = t_init;
while (t < t_final)

  if (timing_run ~= 1)

    % compute exact solution and err
    if (source_term_type == 1)
      T_exact = 1.0 + 0.5*x ...
              + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*kappa*t) ...
              - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*kappa*t) ...
              + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*kappa*t);
    elseif (source_term_type == 2)
      T_exact = 1.0 + 0.5*x ...
              + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*kappa*t) ...
              - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*kappa*t) ...
              + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*kappa*t) ...
              + 50*4/9/pi^2*sin(3/2*pi*x) ...
               *(1-exp(-9/4*pi^2*kappa*t)) ...
              - 70*4/225/pi^2*sin(15/2*pi*x) ...
               *(1-exp(-225/4*pi^2*kappa*t)) ...
              + 100*4/21^2/pi^2*sin(21/2*pi*x) ...
               *(1-exp(-21^2/4*pi^2*kappa*t));
    else
      error('Invalid source term type.');
    end
    err = T-T_exact;
    err(idx_outside) = 0;
    err_L_inf = norm(err,'inf')
    err_L2    = norm(err,2)*dx

    % plot current solution
    figure(1); clf;
    plot(x,T,'bo');
    hold on;
    T_exact_plot = T_exact;  T_exact_plot(idx_outside) = 0; 
    plot(x,T_exact_plot,'r');
    title_string = sprintf('t = %f',t);
    title(title_string);

    % plot current error
    figure(2); clf;
    plot(x,err);
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow

  end

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t; 
  end

  % impose boundary conditions
  if (source_term_type == 1)
    T_bdry = 1.0 + 0.5*x_bdry ...
           + 2*sin(5/2*pi*x_bdry)*exp(-25/4*pi^2*kappa*t) ...
             - 4*sin(11/2*pi*x_bdry)*exp(-121/4*pi^2*kappa*t) ...
           + 3*sin(7/2*pi*x_bdry)*exp(-49/4*pi^2*kappa*t);
  elseif (source_term_type == 2)
    T_bdry = 1.0 + 0.5*x_bdry ...
           + 2*sin(5/2*pi*x_bdry)*exp(-25/4*pi^2*kappa*t) ...
           - 4*sin(11/2*pi*x_bdry)*exp(-121/4*pi^2*kappa*t) ...
           + 3*sin(7/2*pi*x_bdry)*exp(-49/4*pi^2*kappa*t) ...
           + 50*4/9/pi^2*sin(3/2*pi*x_bdry) ...
            *(1-exp(-9/4*pi^2*kappa*t)) ...
           - 70*4/225/pi^2*sin(15/2*pi*x_bdry) ...
            *(1-exp(-225/4*pi^2*kappa*t)) ...
           + 100*4/21^2/pi^2*sin(21/2*pi*x_bdry) ...
            *(1-exp(-21^2/4*pi^2*kappa*t));
  else
    error('Invalid source term type.');
  end

  if (bc_polynomial_interpolant_order == 0)
 
    T(idx_ghostcells) = coefs_bc_bdry_pts.*T_bdry;
 
  elseif (bc_polynomial_interpolant_order == 1)
 
    T(idx_ghostcells) = coefs_bc_bdry_pts.*T_bdry ...
                      + coefs_bc_interp_pts.*T(idx_bc_interp_pts);

  elseif (bc_polynomial_interpolant_order == 2)

    T(idx_ghostcells) = coefs_bc_bdry_pts.*T_bdry ...
                      + coefs_bc_interp_pts(:,1).*T(idx_bc_interp_pts(:,1)) ...
                      + coefs_bc_interp_pts(:,2).*T(idx_bc_interp_pts(:,2));

  elseif (bc_polynomial_interpolant_order == 3)

    T(idx_ghostcells) = coefs_bc_bdry_pts.*T_bdry ...
                      + coefs_bc_interp_pts(:,1).*T(idx_bc_interp_pts(:,1)) ...
                      + coefs_bc_interp_pts(:,2).*T(idx_bc_interp_pts(:,2)) ...
                      + coefs_bc_interp_pts(:,3).*T(idx_bc_interp_pts(:,3));

  else
 
    error('Invalid order for polynomial interpolant.');

  end

  % update solution
  if (source_term_type == 1)
    T = T + dt*kappa*(L*T);
  elseif (source_term_type == 2)
    T = T + dt*kappa*(L*T) + dt*f + 0.5*kappa*dt^2*L*f; 
  else
    error('Invalid source term type.');
  end
  T(idx_outside) = 0;

  % update time
  t = t + dt;

end

% compute exact solution and err
if (source_term_type == 1)
  T_exact = 1.0 + 0.5*x ...
          + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*kappa*t) ...
          - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*kappa*t) ...
          + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*kappa*t);
elseif (source_term_type == 2)
  T_exact = 1.0 + 0.5*x ...
          + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*kappa*t) ...
          - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*kappa*t) ...
          + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*kappa*t) ...
          + 50*4/9/pi^2*sin(3/2*pi*x)*(1-exp(-9/4*pi^2*kappa*t)) ...
          - 70*4/225/pi^2*sin(15/2*pi*x)*(1-exp(-225/4*pi^2*kappa*t)) ...
          + 100*4/21^2/pi^2*sin(21/2*pi*x)*(1-exp(-21^2/4*pi^2*kappa*t));
else
  error('Invalid source term type.');
end
T_exact(idx_outside) = 0;
err = T-T_exact;
err_L_inf = norm(err,'inf')
err_L2    = norm(err,2)*dx

% plot final result
figure(1); clf;
plot(x,T,'bo')
hold on;
plot(x,T_exact,'r')
title_string = sprintf('t = %f',t);
title(title_string);

% plot error
figure(2); clf;
plot(x,err);
title_string = sprintf('t = %f',t);
title(title_string);
