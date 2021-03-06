%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveWaveEqn1dKPY() computes the solutions of the 1d wave equation 
%
%   u_tt = c^2 u_xx + f(x,t)
%
% on the domain -1 < x < 1 subject to the initial conditions
%
%   u(x,0)   = 2*exp( -10 (sin(0.5*pi*x))^2 )
%   u_t(x,0) = -0.5*pi*sin(pi*x) + g(x)
%
% and periodic boundary conditions.
%
% When use_source_term is set to 0, f(x,t) = 0, g(x) = 0 and the analytical 
% solution is given by
% 
%   u = exp(-10*(sin(0.5*pi*(x-c*t))).^2)+exp(-10*(sin(0.5*pi*(x+c*t))).^2)) ...
%     + 0.25/c*(cos(pi*(x+c*t)) - cos(pi*(x-c*t))) 
%
% When use_source_term is set to 1, f(x,t) = 0.125/c*pi^2*sin(pi*(x - 0.25*t)),
% g(x) = 0.5*pi/c*cos(pi*x) and the analytical solution to this problem is 
%
%   u = exp(-10*(sin(0.5*pi*(x-c*t))).^2)+exp(-10*(sin(0.5*pi*(x+c*t))).^2)) ...
%     + 0.25/c*(cos(pi*(x+c*t)) - cos(pi*(x-c*t))) ...
%     - 0.25/c/(c+0.25)*( sin(pi*(x-0.25*t)) - sin(pi*(x+c*t)) ) ...
%     + 0.25/c/(c-0.25)*( sin(pi*(x-0.25*t)) - sin(pi*(x-c*t)) )
%
% The numerical solution is computed on a node-centered grid using the
% direct completely centered discretization of the second-order wave 
% equation introduced by Kreiss, Petersson, and Ystrom (2002).
%
% USAGE:
%   function [u, u_exact, x] = solveWaveEqn1dKPY(c, ...
%                                                use_source_term, ...
%                                                N, dt, ...
%                                                t_final, ...
%                                                debug_on)
%
% Arguments:
% - c:                   wave speed
% - use_source_term:     boolean indicating whether or not the source
%                        term should be included in the wave equation.
%                        The source term is deactivated if use_source_term
%                        is set to 0; otherwise the source term is activated.
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
% CHANGE LOG:
% -----------
% 2008/08:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x] = solveWaveEqn1dKPY(c, ...
                                             use_source_term, ...
                                             N, dt, ...
                                             t_final, ...
                                             debug_on)


% check arguments
if (nargin < 5)
  error('solveWaveEqn1dKPY: missing arguments');
end
if (nargin < 6)
  debug_on = 0;
end

% construct grid
x_lo = -1.0;
x_hi =  1.0;
dx = (x_hi-x_lo)/N;
x = x_lo:dx:x_hi-dx; x = x';  % periodic grid point not included

% construct Laplacian operator (with periodic boundary conditions) 
e = ones(N,1);
L = 1/dx^2*spdiags([e -2*e e], -1:1, N, N);
L(1,N) = 1/dx^2;    % periodic BC at x = -1
L(end,1) = 1/dx^2;  % periodic BC at x =  1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KPY time integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize t
t = 0.0;

% set initial conditions
if (use_source_term > 0)
  u   = 2*exp(-10*(sin(0.5*pi*x)).^2);
  u_t = -0.5*(pi*sin(pi*x) - pi/c*cos(pi*x));
else
  u   = 2*exp(-10*(sin(0.5*pi*x)).^2);
  u_t = -0.5*pi*sin(pi*x);
end

% use second-order forward Euler for first time step
if (use_source_term > 0)
  f = 0.125/c*pi^2*sin(pi*x);
else
  f = 0;
end
u_tt = c^2*L*u + f;
u_next = u + dt*u_t + 0.5*dt^2*u_tt;

% update u_prev and u
u_prev = u;
u = u_next;

% update time
t = t + dt;

while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    if (use_source_term > 0)
      u_exact = (exp(-10*(sin(0.5*pi*(x-c*t))).^2) ...
                +exp(-10*(sin(0.5*pi*(x+c*t))).^2)) ...
              + 0.25/c*(cos(pi*(x+c*t)) - cos(pi*(x-c*t))) ...
              - 0.25/c/(c+0.25)*( sin(pi*(x-0.25*t)) - sin(pi*(x+c*t)) ) ...
              + 0.25/c/(c-0.25)*( sin(pi*(x-0.25*t)) - sin(pi*(x-c*t)) );
    else
      u_exact = (exp(-10*(sin(0.5*pi*(x-c*t))).^2) ...
                +exp(-10*(sin(0.5*pi*(x+c*t))).^2)) ...
              + 0.25/c*(cos(pi*(x+c*t)) - cos(pi*(x-c*t)));
    end
    err = u-u_exact;
    err_L_inf = norm(err,'inf')

    % plot current solution
    figure(1); clf;
    plot(x,u,'bo');
    hold on;
    plot(x,u_exact,'r');
    title_string = sprintf('t = %f',t);
    axis([x_lo x_hi -5 5]);
    title(title_string);

    % plot current error
    figure(2); clf;
    plot(x,err);
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow
    pause

  end %  end case: (debug_on == 1)

  % compute source term
  if (use_source_term > 0)
    f = 0.125/c*pi^2*sin(pi*(x-0.25*t));
  else
    f = 0;
  end

  % compute u_tt
  u_tt = c^2*L*u + f;

  % update solution
  u_next = 2*u - u_prev + dt^2*u_tt;

  % update u_prev and u
  u_prev = u;
  u = u_next;

  % update time
  t = t + dt;

  % if we have overstepped t_final, use linear interpolation to obtain
  % u(t_final).
  % NOTE:  solution is still second-order accurate because error in 
  %        linear interpolation is O(dt^2).
  if (t > t_final)
    dt_final = t-t_final;
    u = dt_final/dt*u_prev + (1-dt_final/dt)*u;
    t = t_final;
  end

end

% compute exact solution 
if (use_source_term > 0)
  u_exact = (exp(-10*(sin(0.5*pi*(x-c*t))).^2) ...
            +exp(-10*(sin(0.5*pi*(x+c*t))).^2)) ...
          + 0.25/c*(cos(pi*(x+c*t)) - cos(pi*(x-c*t))) ...
          - 0.25/c/(c+0.25)*( sin(pi*(x-0.25*t)) - sin(pi*(x+c*t)) ) ...
          + 0.25/c/(c-0.25)*( sin(pi*(x-0.25*t)) - sin(pi*(x-c*t)) );
else
  u_exact = (exp(-10*(sin(0.5*pi*(x-c*t))).^2) ...
            +exp(-10*(sin(0.5*pi*(x+c*t))).^2)) ...
          + 0.25/c*(cos(pi*(x+c*t)) - cos(pi*(x-c*t)));
end

% add back periodic grid point
x = [x; x_hi];
u = [u; u(1)];
u_exact = [u_exact; u_exact(1)];
