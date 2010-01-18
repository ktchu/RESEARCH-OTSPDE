%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqnCrankNicholson1d() computes the solutions of the 1d 
% diffusion equation 
%
%   u_t = D u_xx + f(x)
%
% on the domain 0 < x < 1 with f(x) = 2*t - 42*x^5 - 6*x and boundary
% conditions given by u(0) = 0 and u'(1) = 10.  The initial condition is
% taken to be u(x,0) = x^3 + x^3.  The numerical solution is computed on 
% a node-centered grid using Crank-Nicholson time integration with a 
% second-order central difference approximation to the Laplacian.  A 
% suboptimal time step is used.
%
% USAGE:
%   function [u, u_exact, x, timing_data] = ...
%     solveDiffusionEqnCrankNicholson1d(D, ...
%                                       dx, dt, ...
%                                       t_final, ...
%                                       debug_on, timing_on)
%
% Arguments:
% - D:                   diffusion coefficient
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
% - timing_data:         time used to compute the solution (set up time
%                        excluded).  If timing is not activated, timing_data 
%                        is set to -1.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% CHANGE LOG:
% -----------
% 2008/04:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x, timing_data] = ...
  solveDiffusionEqnCrankNicholson1d(D, ...
                                    dx, dt, ...
                                    t_final, ...
                                    debug_on, timing_on)


% check arguments
if (nargin < 4)
  error('solveDiffusionEqnCrankNicholson1d: missing arguments');
end
if (nargin < 5)
  debug_on = 0;
end
if (nargin < 6)
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
u = sinh(sin(4*pi*x)+cos(3*pi*x)-1);

% Crank-Nicholson time integration
t = 0;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = cos(5*t)*sinh(sin(4*pi*x)+cos(3*pi*x)-1);
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

  % set boundary conditions
  u_0 = 0;
  dudx_1 = 4*pi*cos(5*t)*cosh(2);
  dudx_1 = 0.5*(4*pi*cos(5*t)*cosh(2) + 4*pi*cos(5*(t+dt))*cosh(2));

  % set source term
  f_cur = -5*sin(5*t)*sinh(sin(4*pi*x)+cos(3*pi*x)-1) ...
        - D*cos(5*t)*(4*pi*cos(4*pi*x)-3*pi*sin(3*pi*x)).^2 ...
                    .*sinh(sin(4*pi*x)+cos(3*pi*x)-1) ...
        - D*cos(5*t)*(-16*pi^2*sin(4*pi*x)-9*pi^2*cos(3*pi*x)) ...
                    .*cosh(sin(4*pi*x)+cos(3*pi*x)-1);
  f_next = -5*sin(5*(t+dt))*sinh(sin(4*pi*x)+cos(3*pi*x)-1) ...
         - D*cos(5*(t+dt))*(4*pi*cos(4*pi*x)-3*pi*sin(3*pi*x)).^2 ...
                           .*sinh(sin(4*pi*x)+cos(3*pi*x)-1) ...
         - D*cos(5*(t+dt))*(-16*pi^2*sin(4*pi*x)-9*pi^2*cos(3*pi*x)) ...
                           .*cosh(sin(4*pi*x)+cos(3*pi*x)-1);

  % compute RHS 
  rhs = L_rhs*u + 0.5*dt*(f_cur+f_next);
  rhs(1) = rhs(1) + D*dt/dx/dx*u_0;        % Dirichlet BC at x = 0
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
u_exact = cos(5*t)*sinh(sin(4*pi*x)+cos(3*pi*x)-1);
u_0 = 0;

% add back Dirichlet boundary condition
u = [u_0; u];
u_exact = [u_0; u_exact];
x = [0; x];
