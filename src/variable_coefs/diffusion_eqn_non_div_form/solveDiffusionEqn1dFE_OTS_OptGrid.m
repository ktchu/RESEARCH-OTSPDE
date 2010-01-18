%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqn1dFE_OTS_OptGrid() computes the solutions of the 1d 
% diffusion equation with variable diffusivity
%
%   u_t = D*u_xx + f(x,t)
%
% on the domain 0 < x < 1 subject to the initial conditions
%
%   u(x,0)   = 2 + sin(3*pi*x)
%
% and homogeneous Dirichlet boundary conditions.  The diffusivity and 
% source term are given by
%
%   D(x) = (2+sin(2*pi*x))^2
%
%   f(x,t) = -3*pi*cos(3*pi*(x-t))
%          + 9*pi^2*(2+sin(2*pi*x))^2*sin(3*pi*(x-t))
%
% The analytical solution to this problem is 
%
%   u(x,t) = 2+sin(3*pi*(x-t))
%
% The numerical solution is computed on the variable grid size mesh induced 
% by taking a uniform mesh in a transformed domain where the diffusion 
% coefficient is constant.  Forward Euler time integration is used to 
% advance the solution in time.  The optimal time step of dt = (dy/D_bar)^2/6
% is used.  In this expression, dy is the grid spacing of the uniform mesh
% on the transformed domain and D_bar is the harmonic mean of the square-root
% of the diffusion coefficient.  The expression for D_bar is given below.
%
% On the transformed domain 
%
%   y = 1/pi * ( atan( (2*tan(pi*x)+1)/sqrt(3) ) - pi/6 )
%
% the variable coefficient diffusion equation above becomes a diffusion 
% equation that has a constant coefficient leading-order spatial derivative 
% term:
%
%   u_tt = D_bar^2 u_yy + ...
%
% where D_bar has a value of sqrt(3).  Note that in the change of variables,
% the appropriate branch of atan() must be used to avoid a discontinuity in 
% y.
%
%
% USAGE:
%   function [u, u_exact, x] = solveDiffusionEqn1dFE_OTS_OptGrid(N, ...
%                                                                t_final, ...
%                                                                debug_on)
%
% Arguments:
% - N:                   number of grid cells to use for computation.
%                        NOTE:  num grid points = (num grid cells + 1)
% - t_final:             final time
% - debug_on:            flag indicating whether debugging information
%                        should be displayed.  To turn on debugging,
%                        set debug_on to 1.
%                        (default = 0)
%
% Return values:
% - u:                   numerical solution
% - u_exact:             analytical solution
% - x:                   grid points
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x] = solveDiffusionEqn1dFE_OTS_OptGrid(N, ...
                                                             t_final, ...
                                                             debug_on)


% check arguments
if (nargin < 2)
  error('solveDiffusionEqn1dFE_OTS_OptGrid: missing arguments');
end
if (nargin < 3)
  debug_on = 0;
end

% construct grid
y_lo = 0.0;
y_hi = 1.0;
dy = (y_hi-y_lo)/N;
y = y_lo:dy:y_hi; y = y';  

% compute transformation from y to x
x = 1/pi*atan(0.5*(sqrt(3)*tan(pi*y + pi/6)-1));
x(1) = 0;
idx = find(x < 0);
x(idx) = x(idx) + 1;
x_lo = x(1);
x_hi = x(end);

% construct Laplacian operator (with boundary conditions) 
diag_plus  = (2./(x(3:end)-x(1:end-2))) .* (1./(x(3:end)-x(2:end-1)));
diag_minus = (2./(x(3:end)-x(1:end-2))) .* (1./(x(2:end-1)-x(1:end-2)));
diag_center = -(diag_plus + diag_minus);
diag_plus = [0; 0; diag_plus];
diag_minus = [diag_minus; 0; 0];
diag_center = [0; diag_center; 0];
L = spdiags([diag_minus diag_center diag_plus], -1:1, N+1, N+1);
L(1,:) = 0;    % no need to update boundary condition
L(end,:) = 0;  % no need to update boundary condition

% compute diffusivity and its derivatives
D = (2+sin(2*pi*x)).^2;
D_xx = 8*pi^2*(2*cos(2*pi*x).^2 - 2*sin(2*pi*x) - 1);

% compute optimal time step 
D_bar = sqrt(3);
dt = (dy/D_bar)^2/6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE time integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize t
t = 0.0;

% set initial conditions
u = 2 + sin(3*pi*x);

% initialize f_prev to f(-dt)
f_prev = -3*pi*cos(3*pi*(x+dt)) + 9*pi^2*((2+sin(2*pi*x)).^2).*sin(3*pi*(x+dt));

while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = 2+sin(3*pi*(x-t));
    err = u-u_exact;
    err_L_inf = norm(err,'inf')

    % plot current solution
    figure(1); clf;
    plot(x,u,'bo');
    hold on;
    plot(x,u_exact,'r');
    title_string = sprintf('t = %f',t);
    axis([x_lo x_hi 0 5]);
    title(title_string);

    % plot current error
    figure(2); clf;
    plot(x,err);
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow
    %pause

  end %  end case: (debug_on == 1)

  % compute source term
  f = -3*pi*cos(3*pi*(x-t)) + 9*pi^2*((2+sin(2*pi*x)).^2).*sin(3*pi*(x-t));

  % compute u_xx
  u_xx = L*u;

  % compute u_t
  u_t = D.*u_xx + f;

  % compute correction terms
  f_t = (f-f_prev)/dt;
  f_xx = L*f;

  u_corr = 0.5*dt^2*( D.*D_xx.*u_xx + D.*f_xx + f_t );

  % compute solution at next time 
  u_next = u + dt*u_t + u_corr;

  % update time
  t = t + dt;

  % impose boundary conditions
  u_next(1) = 2 + sin(-3*pi*t);
  u_next(end) = 2 + sin(3*pi*(1-t));

  % update solution
  if (t < t_final)
    u = u_next;
    f_prev = f;

  else
    % if we have overstepped t_final, use linear interpolation in time to 
    % obtain u(t_final).
    % NOTE:  solution is limited to fourth-order accuracy because error in 
    %        linear interpolation is O(dt^2).
    dt_final = t-t_final;
    u = dt_final/dt*u+ (1-dt_final/dt)*u_next;
    t = t_final;
  end

end

% compute exact solution 
u_exact = 2+sin(3*pi*(x-t));
