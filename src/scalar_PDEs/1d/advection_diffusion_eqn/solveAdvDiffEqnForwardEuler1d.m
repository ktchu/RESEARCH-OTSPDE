%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveAdvDiffEqnForwardEuler1d() computes the solutions of the 1d 
% advection-diffusion equation 
%
%   u_t = D u_xx + A u_x
%
% on the domain -10 < x < 10 subject to the initial condition
%
%   u(x,0) = exp( -x^2/4 )
%
% and boundary conditions
%
%   u(-10,t) = 1/sqrt(Dt+1) exp( -0.25*(-10+At)^2/(Dt+1) )
%
%   u(10,t)  = 1/sqrt(Dt+1) exp( -0.25*(10+At)^2/(Dt+1) )
%
% which are derived from the analytical solution of the advection-diffusion
% equation on an infinite domain with a Gaussian initial condition.  The 
% numerical solution is computed on a node-centered grid using forward 
% Euler time integration using second-order central difference 
% approximations for both the advection and the diffusion terms.
%
% USAGE:
%   function [u, u_exact, x] = solveAdvDiffEqnForwardEuler1d( ...
%                                D, A, ...
%                                dx, dt, ...
%                                t_final, ...
%                                debug_on)
%
% Arguments:
% - D:                   diffusion coefficient 
% - A:                   flow speed 
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

function [u, u_exact, x] = solveAdvDiffEqnForwardEuler1d( ...
                             D, A, ...
                             dx, dt, ...
                             t_final, ...
                             debug_on)


% check arguments
if (nargin < 5)
  error('solveAdvDiffEqnForwardEuler1d: missing arguments');
end
if (nargin < 6)
  debug_on = 0;
end

% construct grid
N = 20.0/dx;
x = -10:dx:10; x = x';

% construct Laplacian operator (with boundary conditions) 
e = ones(N+1,1);
L = 1/dx^2*spdiags([e -2*e e], -1:1, N+1, N+1);
L(1,:) = 0;    % no need to compute laplacian for grid point at x = -10
L(end,:) = 0;  % no need to compute laplacian for grid point at x = 10

% construct gradient operator (with boundary conditions) 
G = 0.5/dx*spdiags([-e e], [-1,1], N+1, N+1);
G(1,:) = 0;    % no need to compute gradient for grid point at x = -10
G(end,:) = 0;  % no need to compute gradient for grid point at x = 10
                             
% set initial conditions
u = exp(-0.25*x.^2);


% forward Euler time integration
t = 0.0;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = exp(-0.25*(x+A*t).^2/(D*t+1))/sqrt(D*t+1);
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
  u = u + dt*(D*(L*u) + A*(G*u)); 

  % update time
  t = t + dt;

  % update solution at boundaries
  u(1)   = exp( -0.25*(-10+A*t)^2/(D*t+1) ) / sqrt(D*t+1);
  u(end) = exp( -0.25*(10+A*t)^2/(D*t+1) )  / sqrt(D*t+1);

end

% compute exact solution 
u_exact = exp(-0.25*(x+A*t).^2/(D*t+1))/sqrt(D*t+1); 

