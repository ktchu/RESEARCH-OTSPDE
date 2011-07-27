%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveBurgersEqnForwardEuler1d() computes the solutions of the 1d 
% viscous Burgers equation
%
%   u_t + u u_x = nu u_xx
%
% on the domain 0 < x < 10 subject to the initial condition
%
%   u(x,0) = U + sqrt(nu)*gamma*exp(-(x-U)^2/(4*nu))
%                /( sqrt(pi)*(1 + 0.5*gamma*erfc((x-U)/sqrt(4*nu))) )
%
% and boundary conditions
%
%   u(0,t)  = U + sqrt(nu/T)*gamma*exp(-(-U*T)^2/(4*nu*T))
%                /( sqrt(pi)*(1 + 0.5*gamma*erfc((-U*T)/sqrt(4*nu*T))) )
%
%   u(10,t) = U + sqrt(nu/T)*gamma*exp(-(10-U*T)^2/(4*nu*T))
%                /( sqrt(pi)*(1 + 0.5*gamma*erfc((10-U*T)/sqrt(4*nu*T))) )
%
% which are derived from the analytical single hump solution of the viscous
% Burgers equation on an infinite domain (Whitham, p.102):
%
%   u(x,t) = U + sqrt(nu/T)*gamma*exp(-(x-U*T)^2/(4*nu*T))
%               /( sqrt(pi)*(1 + 0.5*gamma*erfc((x-U*T)/sqrt(4*nu*T))) )
%
% Note that in all of the above expressions, T = t+1, gamma = exp(R)-1, and 
% is the effective Reynolds number R.  The numerical solution is computed 
% on a node-centered grid using forward Euler time integration using 
% second-order central difference approximations for both the nonlinear 
% advection and the diffusion terms.  
%
% USAGE:
%   function [u, u_exact, x, timing_data] = ...
%     solveBurgersEqnForwardEuler1d(nu, U, R, ...
%                                   dx, dt, ...
%                                   t_final, ...
%                                   debug_on, ...
%                                   timing_on)
%
% Arguments:
% - nu:                  viscosity
% - U:                   wave speed
% - R:                   effective Reynolds number
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
% - x:                   grid points
% - timing_data:         array of the following form containing timing data
%                        [total solution time].  If timing is not activated,
%                        timing_data is set to [-1].
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGE LOG:
% -----------
% 2008/02:  Initial version of code. 
% 2011/07:  Added support for collecting timing data.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x, timing_data] = ...
  solveBurgersEqnForwardEuler1d(nu, U, R, ...
                                dx, dt, ...
                                t_final, ...
                                debug_on, ...
                                timing_on)


% check arguments
if (nargin < 6)
  error('solveBurgersEqnForwardEuler1d: missing arguments');
end
if (nargin < 7)
  debug_on = 0;
end
if (nargin < 8)
  timing_on = 0;
end

% construct grid
N = 10.0/dx;
x = 0:dx:10; x = x';

% construct Laplacian operator (with boundary conditions) 
e = ones(N+1,1);
L = 1/dx^2*spdiags([e -2*e e], -1:1, N+1, N+1);
L(1,:) = 0;    % no need to compute laplacian for grid point at left endpt
L(end,:) = 0;  % no need to compute laplacian for grid point at right endpt

% construct gradient operator (with boundary conditions) 
G = 0.5/dx*spdiags([-e e], [-1,1], N+1, N+1);
G(1,:) = 0;    % no need to compute gradient for grid point at left endpt
G(end,:) = 0;  % no need to compute gradient for grid point at right endpt
                             
% compute gamma
gamma = exp(R)-1;

% set initial conditions
u = U +  sqrt(nu)*gamma*exp(-(x-U).^2/(4*nu)) ...
       ./( sqrt(pi)*(1 + 0.5*gamma*erfc((x-U)/sqrt(4*nu))) );

% prepare for collecting timing data
if (timing_on == 1)
  t_start = cputime;
end

% forward Euler time integration
t = 0.0;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    T = t+1;
    u_exact = U +  sqrt(nu/T)*gamma*exp(-(x-U*T).^2/(4*nu*T)) ...
                 ./( sqrt(pi)*(1 + 0.5*gamma*erfc((x-U*T)/sqrt(4*nu*T))) );
    err = u-u_exact;
    err_L_inf = norm(err,'inf')

    % plot current solution
    figure(1); clf;
    plot(x,u,'bo');
    hold on;
    plot(x,u_exact,'r');
    title_string = sprintf('t = %f',t);
    title(title_string);

    % plot current error
    figure(2); clf;
    plot(x,err);
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow

  end %  end case: (debug_on == 1)

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t; 
  end

  % update solution
  u = u + dt*(nu*(L*u) - u.*(G*u)); 

  % update time
  t = t + dt;

  % update solution at boundaries
  T = t+1;
  u(1)   =  U +  sqrt(nu/T)*gamma*exp(-(x(1)-U*T).^2/(4*nu*T)) ...
               ./( sqrt(pi)*(1 + 0.5*gamma*erfc((x(1)-U*T)/sqrt(4*nu*T))) );
  u(end) =  U +  sqrt(nu/T)*gamma*exp(-(x(end)-U*T).^2/(4*nu*T)) ...
             ./( sqrt(pi)*(1 + 0.5*gamma*erfc((x(end)-U*T)/sqrt(4*nu*T))) );

end

% measure time to solve Burgers equation
if (timing_on == 1)
  t_solve = cputime - t_start; 
end

% compute exact solution 
T = t+1;
u_exact = U +  sqrt(nu/T)*gamma*exp(-(x-U*T).^2/(4*nu*T)) ...
             ./( sqrt(pi)*(1 + 0.5*gamma*erfc((x-U*T)/sqrt(4*nu*T))) );

% set timing data
if (timing_on == 1)
  timing_data = [t_solve];
else
  timing_data = [-1];
end
