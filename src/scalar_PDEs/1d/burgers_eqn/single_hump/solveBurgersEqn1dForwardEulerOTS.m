%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveBurgersEqnForwardEulerOTS1d() computes the solutions of the 1d 
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
% advection and the diffusion terms.  The optimal time step of dx^2/(6*nu) 
% and correction terms are used.
%
% USAGE:
%   function [u, u_exact, x] = solveBurgersEqnForwardEulerOTS1d( ...
%                                nu, U, R, ...
%                                dx, ...
%                                t_final, ...
%                                debug_on)
%
% Arguments:
% - nu:                  viscosity
% - U:                   wave speed
% - R:                   effective Reynolds number
% - dx:                  grid spacing
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
% 2008/02:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x] = solveBurgersEqnForwardEulerOTS1d( ...
                             nu, U, R, ...
                             dx, ...
                             t_final, ...
                             debug_on)


% check arguments
if (nargin < 5)
  error('solveBurgersEqnForwardEulerOTS1d: missing arguments');
end
if (nargin < 6)
  debug_on = 0;
end

% construct grid
N = 10.0/dx;
x = 0:dx:10; x = x';

% compute optimal time step
dt = dx^2/6/nu;

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
  u = u + dt*(nu*(L*u) - u.*(G*u)) ...
    - 0.5*dt^2*( 4*nu*(G*u).*(L*u) - 2*u.*(G*u).^2 - (u.^2).*(L*u) ); 

  % update time
  t = t + dt;

  % update solution at boundaries
  T = t+1;
  u(1)   =  U +  sqrt(nu/T)*gamma*exp(-(x(1)-U*T).^2/(4*nu*T)) ...
               ./( sqrt(pi)*(1 + 0.5*gamma*erfc((x(1)-U*T)/sqrt(4*nu*T))) );
  u(end) =  U +  sqrt(nu/T)*gamma*exp(-(x(end)-U*T).^2/(4*nu*T)) ...
             ./( sqrt(pi)*(1 + 0.5*gamma*erfc((x(end)-U*T)/sqrt(4*nu*T))) );

end

% compute exact solution 
T = t+1;
u_exact = U +  sqrt(nu/T)*gamma*exp(-(x-U*T).^2/(4*nu*T)) ...
             ./( sqrt(pi)*(1 + 0.5*gamma*erfc((x-U*T)/sqrt(4*nu*T))) );
