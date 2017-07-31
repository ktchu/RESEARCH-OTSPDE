%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveAdvectionEqn1dForwardEulerOTS() computes the solutions of the 
% system of 1d advection equations
%
%   u_t = -3 u_x
%   v_t =  2 v_x
%
% on the domain -10 < x < 10 subject to the initial conditions
%
% u(x,0) = exp( -x^2/4)
% v(x,0) = sin(0.5*pi*x)
%
% and boundary conditions
%
% u(-10, t) = exp( -(-10-3*t)^2/4 )
% u( 10, t) = exp( -(10-3*t)^2/4 )
%
% v(-10, t) = sin(0.5*pi*(-10+2*t))
% v( 10, t) = sin(0.5*pi*(10+2*t))
%
% The numerical solution is computed on a node-centered grid using 
% forward Euler time integration with a first-order upwind difference 
% approximation for the gradients.  Optimal time stepping is used in 
% conjunction with a quadratic time interpolation scheme to boost the
% order of accuracy of the scheme from first- to second-order.
%
% USAGE:
%   function [u, v, u_exact, v_exact, x, timing_data] = ...
%     solveAdvectionEqn1dForwardEulerOTS(dx, ...
%                                        t_init, t_final, ...
%                                        debug_on, timing_on)
%
% Arguments:
% - dx:                  grid spacing
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
% 2008/07/31:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, v, u_exact, v_exact, x, timing_data] = ...
  solveAdvectionEqn1dForwardEulerOTS(dx, ...
                                     t_init, t_final, ...
                                     debug_on, timing_on)


% check arguments
if (nargin < 3)
  error('solveAdvectionEqn1dForwardEulerOTS: missing arguments');
end
if (nargin < 4)
  debug_on = 0;
end
if (nargin < 5)
  timing_on = 0;
end

% construct grid
x_hi = 10;
x_lo = -10;
N = (x_hi-x_lo)/dx + 1;
x = x_lo:dx:x_hi;  x = x';

% compute optimal time steps
dt_u = dx/3;
dt_v = dx/2;
dt = min(dt_u, dt_v);

% construct gradient operator (with boundary conditions) 
e = ones(N,1);                         
G_plus = 1/dx*spdiags([-e e], 0:1, N, N);
G_plus(1,1) = 0; G_plus(1,2) = 0;   % no need to compute gradient at boundary 
G_plus(N,N-1) = 0; G_plus(N,N) = 0; % no need to compute gradient at boundary 
G_minus = 1/dx*spdiags([-e e], -1:0, N, N);
G_minus(1,1) = 0; G_minus(1,2) = 0;   % no need to compute gradient at boundary 
G_minus(N,N-1) = 0; G_minus(N,N) = 0; % no need to compute gradient at boundary 

% start clock to measure computation 
if (timing_on == 1)
  t_start = cputime;
end

% set initial conditions
u = exp(-x.^2/4);
v = sin(0.5*pi*x);

% forward Euler time integration
t = t_init;

% take first time step using suboptimal time step
% NOTE: this only contributes an O(dt^2) error to the final solution.

% save solution at previous step
u_prev = u;
v_prev = v;

% advance each component of solution 
u = u - 3*dt*G_minus*u; 
v = v + 2*dt*G_plus*v; 

% update time
t = t + dt;

% fix boundary conditions
u(1) = exp(-(x_lo-3*t)^2/4);
u(N) = exp(-(x_hi-3*t)^2/4);
v(1) = sin(0.5*pi*(x_lo+2*t));
v(N) = sin(0.5*pi*(x_hi+2*t));

% take remaining time steps - synchronization uses quadratic interpolation 
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = exp(-(x-3*t).^2/4);
    v_exact = sin(0.5*pi*(x+2*t));
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

  if (t+dt < t_final)
    % advance solution using optimal time step and 
    % quadratic interpolation in time

    % impose boundary conditions
    u(1) = exp(-(x_lo-3*t)^2/4);
    u(N) = exp(-(x_hi-3*t)^2/4);
    v(1) = sin(0.5*pi*(x_lo+2*t));
    v(N) = sin(0.5*pi*(x_hi+2*t));

    % advance each component of solution using optimal time steps
    u_ots= u - 3*dt_u*G_minus*u; 
    v_ots= v + 2*dt_v*G_plus*v; 

    % synchronize solutions using quadratic interpolation 
    u_next = (dt-dt_u)/(dt+dt_u)*u_prev - 2*(dt-dt_u)/dt_u*u ...
           + 2*dt^2/dt_u/(dt+dt_u)*u_ots;
    v_next = (dt-dt_v)/(dt+dt_v)*v_prev - 2*(dt-dt_v)/dt_v*v ...
           + 2*dt^2/dt_v/(dt+dt_v)*v_ots;

    % save solution at previous time step and update the current solution
    u_prev = u;
    v_prev = v;
    u = u_next;
    v = v_next;

  else
    % adjust dt so we don't overstep t_final
    dt = t_final - t; 

    % advance solution using suboptimal time step.
    % NOTE: this only contributes an O(dt^2) error to the final solution.
    u = u - 3*dt*G_minus*u; 
    v = v + 2*dt*G_plus*v; 

  end

  % update time
  t = t + dt;

  % fix boundary conditions
  u(1) = exp(-(x_lo-3*t)^2/4);
  u(N) = exp(-(x_hi-3*t)^2/4);
  v(1) = sin(0.5*pi*(x_lo+2*t));
  v(N) = sin(0.5*pi*(x_hi+2*t));

end

% measure time to solve advection equations
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
u_exact = exp(-(x-3*t).^2/4);
v_exact = sin(0.5*pi*(x+2*t));
