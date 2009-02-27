%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqn1dFE_OTS_OptGrid() computes the solutions of the 1d 
% diffusion equation with variable diffusivity
%
%   u_t = d^2/dx^2 ( D(x) * u ) + f(x,t)
%       = D*u_xx + 2*D_x*u_x + D_xx*u + f(x,t)
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
%   f(x,t) = -pi*(- 16*pi 
%                 - 32*pi*sin(2*pi*x)
%                 + 32*pi*cos(2*pi*x)^2
%                 + 3*cos(3*pi*(x-t))
%                 - 53*pi*sin(3*pi*(x-t))
%                 + 48*pi*cos(2*pi*x) * cos(3*pi*(x-t)) 
%                 - 52*pi*sin(2*pi*x) * sin(3*pi*(x-t)) 
%                 + 25*pi*cos(2*pi*x)^2 * sin(3*pi*(x-t))
%                 + 24*pi*cos(2*pi*x)*sin(2*pi*x) * cos(3*pi*(x-t)) )
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

% construct gradient operator (with boundary conditions) 
diag_plus  =  1./(x(3:end)-x(1:end-2)); 
diag_minus = -diag_plus;
diag_plus = [0; 0; diag_plus];
diag_minus = [diag_minus; 0; 0];
G = spdiags([diag_minus diag_plus], [-1,1], N+1, N+1);
G(1,:) = 0;    % no need to update boundary condition
G(end,:) = 0;  % no need to update boundary condition

% compute diffusivity and its derivatives
D = (2+sin(2*pi*x)).^2;
D_x = 4*pi*(2 + sin(2*pi*x)).*cos(2*pi*x);
D_xx = 8*pi^2*(2*cos(2*pi*x).^2 - 2*sin(2*pi*x) - 1);
D_xxx = -32*pi^3*cos(2*pi*x).*(2*sin(2*pi*x) + 1);
D_xxxx = 64*pi^4*(2 - 4*cos(2*pi*x).^2 + sin(2*pi*x));

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

% initialize f_prev
f_prev = -pi*(- 16*pi ...
              - 32*pi*sin(2*pi*x) ...
              + 32*pi*cos(2*pi*x).^2 ...
              + 3*cos(3*pi*(x+dt)) ...
              - 53*pi*sin(3*pi*(x+dt)) ...
              + 48*pi*cos(2*pi*x) .* cos(3*pi*(x+dt))  ...
              - 52*pi*sin(2*pi*x) .* sin(3*pi*(x+dt))  ...
              + 25*pi*cos(2*pi*x).^2 .* sin(3*pi*(x+dt)) ...
              + 24*pi*cos(2*pi*x).*sin(2*pi*x) .* cos(3*pi*(x+dt)) );

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
  f = -pi*(- 16*pi ...
           - 32*pi*sin(2*pi*x) ...
           + 32*pi*cos(2*pi*x).^2 ...
           + 3*cos(3*pi*(x-t)) ...
           - 53*pi*sin(3*pi*(x-t)) ...
           + 48*pi*cos(2*pi*x) .* cos(3*pi*(x-t))  ...
           - 52*pi*sin(2*pi*x) .* sin(3*pi*(x-t))  ...
           + 25*pi*cos(2*pi*x).^2 .* sin(3*pi*(x-t)) ...
           + 24*pi*cos(2*pi*x).*sin(2*pi*x) .* cos(3*pi*(x-t)) );

%  f_t = pi^2*(- 9*sin(3*pi*(x-t)) ...
%              - 159*pi*cos(3*pi*(x-t)) ...
%              - 144*pi*cos(2*pi*x) .* sin(3*pi*(x-t))  ...
%              - 156*pi*sin(2*pi*x) .* cos(3*pi*(x-t))  ...
%              + 75*pi*cos(2*pi*x).^2 .* cos(3*pi*(x-t)) ...
%              - 72*pi*cos(2*pi*x).*sin(2*pi*x) .* sin(3*pi*(x-t)) );

%  f_x = -pi^2*(- 64*pi*cos(2*pi*x) ...
%               - 64*pi*sin(4*pi*x) ...
%               - 9*sin(3*pi*(x-t)) ...
%               - 207*pi*cos(3*pi*(x-t)) ...
%               - 252*pi*sin(2*pi*x) .* cos(3*pi*(x-t))  ...
%               - 248*pi*cos(2*pi*x) .* sin(3*pi*(x-t))  ...
%               + 171*pi*cos(2*pi*x).^2 .* cos(3*pi*(x-t)) ...
%               - 86*pi*sin(4*pi*x) .* sin(3*pi*(x-t)) );

%  f_xx = pi^3*(-256*pi ...
%               + 27*cos(3*pi*(x-t)) ...
%               + 512*pi*cos(2*pi*x).^2 ...
%               - 965*pi*sin(3*pi*(x-t)) ...
%               +1201*pi*cos(2*pi*x).^2.*sin(3*pi*(x-t)) ...
%               +1200*pi*cos(2*pi*x).*sin(2*pi*x).*cos(3*pi*(x-t)) ...
%               +1248*pi*cos(2*pi*x) .* cos(3*pi*(x-t)) ...
%               -1252*pi*sin(2*pi*x) .* sin(3*pi*(x-t)) ...
%               - 128*pi*sin(2*pi*x) );

  % compute u_x and u_xx
  u_x  = G*u;
  u_xx = L*u;

  % compute u_t
  u_t = D.*u_xx + 2*D_x.*u_x + D_xx.*u + f;

  % compute correction terms
  f_t = (f-f_prev)/dt;
  f_x = G*f;
  f_xx = L*f;

  u_corr = 0.5*dt^2 ...
         * ( D.*D_xxxx.*u + 4*D.*D_xxx.*u_x + 6*D.*D_xx.*u_xx + D.*f_xx ...
           + 2*D_x.*D_xxx.*u + 6*D_x.*D_xx.*u_x + 2*D_x.*f_x ...
           + D_xx.*u_t + f_t );

  % compute solution at next time 
  u_next = u + dt*u_t + u_corr;

  % update f_prev
  f_prev = f;

  % update time
  t = t + dt;

  % impose boundary conditions
  u_next(1) = 2 + sin(-3*pi*t);
  u_next(end) = 2 + sin(3*pi*(1-t));

  % update solution
  if (t < t_final)
    u = u_next;

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
