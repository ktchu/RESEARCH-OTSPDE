%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqn1dFE() computes the solutions of the 1d diffusion 
% equation with variable diffusivity
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
% The numerical solution is computed on a node-centered grid using 
% central difference operators for the Laplacian and gradient term of
% the expanded equation and forward Euler time integration.
%
% USAGE:
%   function [u, u_exact, x] = solveDiffusionEqn1dFE(N, dt, ...
%                                                    t_final, ...
%                                                    debug_on)
%
% Arguments:
% - N:                   number of grid cells to use for computation.
%                        NOTE:  num grid points = (num grid cells + 1)
% - dt:                  time step
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

function [u, u_exact, x] = solveDiffusionEqn1dFE(N, dt, ...
                                                 t_final, ...
                                                 debug_on)


% check arguments
if (nargin < 3)
  error('solveDiffusionEqn1dFE: missing arguments');
end
if (nargin < 4)
  debug_on = 0;
end

% construct grid
x_lo = 0.0;
x_hi = 1.0;
dx = (x_hi-x_lo)/N;
x = x_lo:dx:x_hi; x = x';  

% construct Laplacian operator (with boundary conditions) 
e = ones(N+1,1);
L = 1/dx^2*spdiags([e -2*e e], -1:1, N+1, N+1);
L(1,:) = 0;
L(end,:) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE time integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize t
t = 0.0;

% set initial conditions
u = 2 + sin(3*pi*x);

% compute diffusivity and its derivative
D = (2+sin(2*pi*x)).^2;
D_x = 4*pi*(2 + sin(2*pi*x)).*cos(2*pi*x);
D_xx = 8*pi^2*(2*cos(2*pi*x).^2 - 2*sin(2*pi*x) - 1);


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

  % compute u_t
  u_t = D.*(L*u) + f;

  % compute solution at next time 
  u_next = u + dt*u_t;

  % impose boundary conditions
  u_next(1) = 2 + sin(-3*pi*t);
  u_next(end) = 2 + sin(3*pi*(1-t));

  % update time
  t = t + dt;

  % update solution
  if (t < t_final)
    u = u_next;

  else
    % if we have overstepped t_final, use linear interpolation to obtain
    % u(t_final).
    dt_final = t-t_final;
    u = dt_final/dt*u + (1-dt_final/dt)*u_next;
    t = t_final;

  end

end

% compute exact solution 
u_exact = 2+sin(3*pi*(x-t));
