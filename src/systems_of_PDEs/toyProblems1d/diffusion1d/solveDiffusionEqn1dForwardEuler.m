%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqn1dForwardEuler() computes the solutions of the system 
% of 1d diffusion equations
%
%   u_t = (1/7) u_xx
%   v_t = (1/4) v_xx
%
% on the domain -1 < x < 1 subject to the initial conditions
%
% u(x,0) = 1 + sin(4*pi*x)
% v(x,0) = 3*cos(2*pi*x) + sin(3*pi*x)
%
% and boundary conditions
%
% u(-1, t) = 1
% u( 1, t) = 1
%
% v(-1, t) = 3*exp(-pi^2*t)
% v( 1, t) = 3*exp(-pi^2*t)
%
% The analytical solution to this system of equations is given by
%
% u(x,t) = 1 + exp(-16/7*pi^2*t)*sin(4*pi*x)
% v(x,t) = 3*exp(-pi^2*t)*cos(2*pi*x) + exp(-9/4*pi^2*t)*sin(3*pi*x)
%
% The numerical solution is computed on a node-centered grid using 
% forward Euler time integration with a first-order upwind difference 
% approximation for the Laplacians.
%
% USAGE:
%   function [u, v, u_exact, v_exact, x, timing_data] = ...
%     solveDiffusionEqn1dForwardEuler(dx, dt, ...
%                                     t_init, t_final, ...
%                                     debug_on, timing_on)
%
% Arguments:
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
% - u:                   numerical solution for u
% - v:                   numerical solution for v
% - u_exact:             analytical solution for u
% - v_exact:             analytical solution for v
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
% 2008/08/03:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, v, u_exact, v_exact, x, timing_data] = ...
  solveDiffusionEqn1dForwardEuler(dx, dt, ...
                                  t_init, t_final, ...
                                  debug_on, timing_on)


% check arguments
if (nargin < 4)
  error('solveDiffusionEqn1dForwardEuler: missing arguments');
end
if (nargin < 5)
  debug_on = 0;
end
if (nargin < 6)
  timing_on = 0;
end

% construct grid
x_hi = 1;
x_lo = -1;
N = (x_hi-x_lo)/dx + 1;
x = x_lo:dx:x_hi;  x = x';

% construct Laplacian operator (with boundary conditions) 
e = ones(N,1);                         
L = 1/dx^2*spdiags([e -2*e e], -1:1, N, N);
L(1,1) = 0; L(1,2) = 0;   % no need to compute Laplacian at x = -1
L(N,N-1) = 0; L(N,N) = 0; % no need to compute Laplacian at x = 1

% start clock to measure computation 
if (timing_on == 1)
  t_start = cputime;
end

% set initial conditions
u = 1 + sin(4*pi*x);
v = 3*cos(2*pi*x) + sin(3*pi*x);

% forward Euler time integration
t = t_init;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = 1 + exp(-16/7*pi^2*t)*sin(4*pi*x);
    v_exact = 3*exp(-pi^2*t)*cos(2*pi*x) + exp(-9/4*pi^2*t)*sin(3*pi*x);
    err_u = u-u_exact;
    err_u_L_inf = norm(err_u,'inf')
    err_v = v-v_exact;
    err_v_L_inf = norm(err_v,'inf')

    % plot current solution
    figure(1); clf;
    plot(x,u,'bo');
    hold on;
    plot(x,u_exact,'r');
    plot(x,v,'cs');
    plot(x,v_exact,'g');
    title_string = sprintf('t = %f',t);
    title(title_string);

    % plot current error
    figure(2); clf;
    plot(x,err_u,'b');
    hold on;
    plot(x,err_v,'g');
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow

  end %  end case: (debug_on == 1)

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t; 
  end

  % impose boundary conditions
  u(1) = 1;
  u(N) = 1;
  v(1) = 3*exp(-pi^2*t);
  v(N) = 3*exp(-pi^2*t);

  % update solution
  u = u + dt*L*u/7; 
  v = v + dt*L*v/4; 

  % update time
  t = t + dt;

  % fix boundary conditions
  u(1) = 1;
  u(N) = 1;
  v(1) = 3*exp(-pi^2*t);
  v(N) = 3*exp(-pi^2*t);

end

% measure time to solve diffusion equations
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
    comp_t_solve = sprintf('  Solution Time: %f', t_solve);
    disp(comp_time_str);
    disp('==================================================================');
  end

else

  timing_data = -1;

end

% compute exact solution 
u_exact = 1 + exp(-16/7*pi^2*t)*sin(4*pi*x);
v_exact = 3*exp(-pi^2*t)*cos(2*pi*x) + exp(-9/4*pi^2*t)*sin(3*pi*x);
