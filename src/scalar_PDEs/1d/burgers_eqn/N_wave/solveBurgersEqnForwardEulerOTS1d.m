%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveBurgersEqnForwardEulerOTS1d() computes the solutions of the 1d 
% viscous Burgers equation
%
%   u_t + u u_x = nu u_xx
%
% on the domain -10 < x < 10 subject to the initial condition
%
%   u(x,0) = sqrt(nu) * (exp(R)-1) * exp(-(x-U)^2/(4*nu))
%          / ( sqrt(pi) + (exp(R)-1) * erfc((x-U)/sqrt(4*nu)) ) + U
%
% and boundary conditions
%
%   u(-10,t)  = sqrt(nu/T) * (exp(R)-1) * exp(-(-10-U*T)^2/(4*nu*T))
%             / ( sqrt(pi) + (exp(R)-1) * erfc(-10-U*T/sqrt(4*nu*T)) ) + U
%
%   u(10,t)  = sqrt(nu/T) * (exp(R)-1) * exp(-(10-U*T)^2/(4*nu*T))
%            / ( sqrt(pi) + (exp(R)-1) * erfc((10-U*T)/sqrt(4*nu*T)) ) + U
%
% which are derived from the analytical N-wave solution of the viscous
% Burgers equation on an infinite domain (Whitham, p.107):
%
%   u(x,t) = sqrt(nu/T) * (exp(R)-1) * exp(-(x-U*T)^2/(4*nu*T))
%          / ( sqrt(pi) + (exp(R)-1) * erfc((x-U*T)/sqrt(4*nu*T)) ) + U
%
% Note that in all of the above expressions, T = t+1 and the Reynolds number
% R is set to 2.0/nu.  The numerical solution is computed on a node-centered 
% grid using forward Euler time integration using second-order central 
% difference approximations for both the nonlinear advection and the diffusion 
% terms.  The optimal time step of dx^2/(6*nu) and correction terms are used.
%
% USAGE:
%   function [u, u_exact, x] = solveBurgersEqnForwardEulerOTS1d( ...
%                                nu, ...
%                                dx, ...
%                                t_final, ...
%                                debug_on)
%
% Arguments:
% - nu:                  viscosity
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
                             nu, ...
                             dx, ...
                             t_final, ...
                             debug_on)


% check arguments
if (nargin < 3)
  error('solveBurgersEqnForwardEulerOTS1d: missing arguments');
end
if (nargin < 4)
  debug_on = 0;
end

% construct grid
N = 20.0/dx;
x = -10:dx:10; x = x';

% compute optimal time step
dt = dx^2/6/nu;

% construct Laplacian operator (with boundary conditions) 
L = (diag(ones(N,1),1) - 2*diag(ones(N+1,1),0) + diag(ones(N,1),-1))/dx/dx;
L(1,:) = 0;    % no need to compute laplacian for grid point at left endpt
L(end,:) = 0;  % no need to compute laplacian for grid point at right endpt

% construct gradient operator (with boundary conditions) 
G = (diag(ones(N,1),1) - diag(ones(N,1),-1))/(2*dx);
G(1,:) = 0;    % no need to compute gradient for grid point at left endpt
G(end,:) = 0;  % no need to compute gradient for grid point at right endpt

% compute Reynolds number
R = 2.0/nu; 
R = 0.5/nu; 
                             
% set initial conditions
u =  x .* exp(-x.^2/(4*nu)) ./ (1 + exp(-x.^2/(4*nu)));

% forward Euler time integration
t = 0.0;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    T = t+1;
    u_exact  =  x/T .* (sqrt(1/T)*exp(-x.^2/(4*nu*T))) ...
                    ./ (1 + sqrt(1/T)*exp(-x.^2/(4*nu*T)));
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

  % update soluation at boundaries
  T = t+1;
  u(1)   =  x(1)/T * sqrt(1/T)*exp(-x(1)^2/(4*nu*T)) ...
                   / (1 + sqrt(1/T)*exp(-x(1)^2/(4*nu*T)));
  u(end) =  x(end)/T * sqrt(1/T)*exp(-x(end)^2/(4*nu*T)) ...
                     / (1 + sqrt(1/T)*exp(-x(end)^2/(4*nu*T)));


end

% compute exact solution 
T = t+1;
u_exact  =  x/T .* (sqrt(1/T)*exp(-x.^2/(4*nu*T))) ...
                ./ (1 + sqrt(1/T)*exp(-x.^2/(4*nu*T)));
