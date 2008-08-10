%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqnCrankNicholson1d() computes the solutions of the 1d
% diffusion equation
%
%   u_t = D u_xx + f(x)
%
% on the domain 0 < x < 1 with a Dirichlet boundary condition imposed at
% x = 0 and a Neumann boundary condition imposed at x = 1.  The numerical
% solution is computed on a node-centered grid using Crank-Nicholson time
% integration using the standard second-order discretization of the 
% Laplacian.  A suboptimal time step is used.
%
% USAGE:
%   function [u, u_exact, x, timing_data] = ...
%     solveDiffusionEqnCrankNicholson1d(D, ...
%                                       source_term_type, ...
%                                       u_0, dudx_1, ...
%                                       dx, dt, ...
%                                       t_init, t_final, ...
%                                       debug_on, timing_on)
%
% Arguments:
% - D:                   diffusion coefficient
% - source_term_type:    type of source term to use in computation.
%                        0:  f = 0
%                        1:
%
%           f = 5*sin(3/2*pi*x) - 7*sin(15/2*pi*x) + 10*sin(21/2*pi*x)
%
% - u_0:                 Dirichlet boundary condition at x = 0
% - dudx_1:              Neumann boundary condition at x = 1
% - dx:                  grid spacing
% - dt:                  time step
% - t_init:              initial time
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
% - timing_data:         time used to compute the solution (set up time
%                        excluded).  If timing is not activated, timing_data 
%                        is set to -1.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% CHANGE LOG:
% -----------
% 2007/09:  Initial version of code. 
% 2008/02:  Changed code to be a function.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x, timing_data] = ...
  solveDiffusionEqnCrankNicholson1d(D, ...
                                    source_term_type, ...
                                    u_0, dudx_1, ...
                                    dx, dt, ...
                                    t_init, t_final, ...
                                    debug_on, timing_on)


% check arguments
if (nargin < 8)
  error('solveDiffusionEqnCrankNicholson1d: missing arguments');
end
if (nargin < 9)
  debug_on = 0;
end
if (nargin < 10)
  timing_on = 0;
end

% construct grid
N = 1/dx + 1;
x = dx:dx:1;  x = x';   % exclude Dirichlet boundary
x_plot = [0; x];

% construct Laplacian operators (with boundary conditions) 
e = ones(N,1);                         
L = 1/dx^2*spdiags([e -2*e e], -1:1, N-1, N-1);
L(N-1,N-2) = 2/dx/dx; L(N-1,N-1) = -2/dx/dx; % modified for Neumann BC at x = 1

L_rhs = speye(N-1) + D*dt/2*L;
L_lhs = speye(N-1) - D*dt/2*L;

% start clock to measure computation 
if (timing_on == 1)
  t_start = cputime;
end

% set initial conditions
u = 1.0 + 0.5*x + 2*sin(5/2*pi*x) - 4*sin(11/2*pi*x) + 3*sin(7/2*pi*x);

% set source term
if (source_term_type == 0)
  f = zeros(size(x));
elseif (source_term_type == 1)
  f = 5*sin(3/2*pi*x) - 7*sin(15/2*pi*x) + 10*sin(21/2*pi*x);
else
  error('solveDiffusionEqnForwardEuler1d: Invalid source term type. Valid valuesare 0 and 1');
end

% Crank-Nicholson time integration
t = t_init;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    if (f == 0)
      u_exact = 1.0 + 0.5*x ...
              + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*D*t) ...
              - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*D*t) ...
              + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*D*t);
    else
      u_exact = 1.0 + 0.5*x ...
              + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*D*t) ...
              - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*D*t) ...
              + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*D*t) ...
              + 5/D*4/9/pi^2*sin(3/2*pi*x)*(1-exp(-9/4*pi^2*D*t)) ...
              - 7/D*4/225/pi^2*sin(15/2*pi*x)*(1-exp(-225/4*pi^2*D*t)) ...
              + 10/D*4/21^2/pi^2*sin(21/2*pi*x)*(1-exp(-21^2/4*pi^2*D*t));
    end
    err = u-u_exact;
    err_L_inf = norm(err,'inf')

    % plot current solution
    figure(1); clf;
    plot(x_plot,[u_0; u],'bo');
    hold on;
    plot(x_plot,[u_0; u_exact],'r');
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
    L_rhs = eye(N-1) + D*dt/2*L;
    L_lhs = eye(N-1) - D*dt/2*L;
  end

  % compute RHS 
  rhs = L_rhs*u + dt*f;
  rhs(1) = rhs(1) + D*dt/dx/dx*u_0;      % Dirichlet BC at x = 0
  rhs(N-1) = rhs(N-1) + 2*D*dt/dx*dudx_1;  % Neumann BC at x = 1

  % update solution
  u = L_lhs\rhs;

  % update time
  t = t + dt;

end

% measure time to solve diffusion equation
if (timing_on == 1)
  t_solve = cputime - t_start;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output timing statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (timing_on)
  timing_data = t_solve;

  if (debug_on)
    disp(' ');
    disp('==================================================================');
    disp('Computation Statistics');
    comp_time_str = sprintf('  Solution Time: %f', t_solve);
    disp(comp_time_str);
    disp('==================================================================');
  end

else

  timing_data = -1;

end

% compute exact solution
if (source_term_type == 0)
  u_exact = 1.0 + 0.5*x ...
          + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*D*t) ...
          - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*D*t) ...
          + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*D*t);
else      
  u_exact = 1.0 + 0.5*x ...
          + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*D*t) ...
          - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*D*t) ...
          + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*D*t) ...
          + 5/D*4/9/pi^2*sin(3/2*pi*x)*(1-exp(-9/4*pi^2*D*t)) ...
          - 7/D*4/225/pi^2*sin(15/2*pi*x)*(1-exp(-225/4*pi^2*D*t)) ...
          + 10/D*4/21^2/pi^2*sin(21/2*pi*x)*(1-exp(-21^2/4*pi^2*D*t));
end

% add back Dirichlet boundary condition
u = [u_0; u];
u_exact = [u_0; u_exact];
x = [0; x];
