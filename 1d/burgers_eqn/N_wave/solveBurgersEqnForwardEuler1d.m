%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveBurgersEqnForwardEuler1d() computes the solutions of the 1d 
% viscous Burgers equation
%
%   u_t + u u_x = nu u_xx
%
% on the domain -10 < x < 10 subject to the initial condition
%
%   u(x,0) =  x * exp(-x^2/(4*nu)) / (1 + exp(-x^2/(4*nu)))
%
% and boundary conditions
%
%   u(-10,t) = -10/T * (sqrt(1/T)*exp(-100/(4*nu*T))) ...
%                    / (1 + sqrt(1/T)*exp(-100^2/(4*nu*T)))
%   u(10,t)  = 10/T * (sqrt(1/T)*exp(-100/(4*nu*T))) ...
%                   / (1 + sqrt(1/T)*exp(-100^2/(4*nu*T)))
%
% which are derived from the analytical N-wave solution of the viscous
% Burgers equation on an infinite domain (Whitham, p.107):
%
%   u(x,t) = x/T * (sqrt(1/T)*exp(-x^2/(4*nu*T))) ...
%                / (1 + sqrt(1/T)*exp(-x^2/(4*nu*T)))
%
% Note that in all of the above expressions, T = t+1.  The numerical solution 
% is computed on a node-centered grid using forward Euler time integration 
% using second-order central difference approximations for both the nonlinear 
% advection and the diffusion terms.
%
% USAGE:
%   function [u, u_exact, x] = solveBurgersEqnForwardEuler1d( ...
%                                nu, ...
%                                dx, dt, ...
%                                t_final, ...
%                                debug_on)
%
% Arguments:
% - nu:                  viscosity
% - dx:                  grid spacing
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
% 2008/02:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x] = solveBurgersEqnForwardEuler1d( ...
                             nu, ...
                             dx, dt, ...
                             t_final, ...
                             debug_on)


% check arguments
if (nargin < 4)
  error('solveBurgersEqnForwardEuler1d: missing arguments');
end
if (nargin < 5)
  debug_on = 0;
end

% construct grid
N = 20.0/dx;
x = -10:dx:10; x = x';

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
  u = u + dt*(nu*(L*u) - u.*(G*u)); 

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
