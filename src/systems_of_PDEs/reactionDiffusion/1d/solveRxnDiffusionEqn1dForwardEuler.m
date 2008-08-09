%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveRxnDiffusionEqn1dForwardEuler() computes the solutions of the system 
% of 1d reaction-diffusion equations
%
%   u_t = u_xx   + 2*(v^3 - 3*v^2*u + v*u^2 + u^3)
%   v_t = D*v_xx + ( (x-t)*v^3 - u^3 + 6*D*v*u^2 - 2*D*v^3 )
%
% on the domain -10 < x < 10 subject to the initial conditions
%
% u(x,0) = 1/(1+x^2)
% v(x,0) = x/(1+x^2)
%
% and boundary conditions
%
% u(-10,t) = 1/(1+(-10-t)^2)
% u( 10,t) = 1/(1+(10-t)^2)
%
% v(-10,t) = (-10-t)/(1+(-10-t)^2)
% v( 10,t) = (10-t)/(1+(10-t)^2)
%
% The analytical solution to this system of equations is given by
%
% u(x,t) = 1/(1+(x-t)^2)
% v(x,t) = (x-t)/(1+(x-t)^2)
% 
% The numerical solution is computed on a node-centered grid using 
% forward Euler time integration with a first-order upwind difference 
% approximation for the Laplacians.
%
% USAGE:
%   function [u, v, u_exact, v_exact, x, timing_data] = ...
%     solveRxnDiffusionEqn1dForwardEuler(D, ...
%                                        N, dt, ...
%                                        t_init, t_final, ...
%                                        debug_on, timing_on)
%
% Arguments:
% - D:                   diffusion constant for v 
% - N:                   number of grid cells to use for computation.
%                        NOTE:  num grid points = (num grid cells + 1)
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
% 2008/08/08:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, v, u_exact, v_exact, x, timing_data] = ...
  solveRxnDiffusionEqn1dForwardEuler(D, ...
                                     N, dt, ...
                                     t_init, t_final, ...
                                     debug_on, timing_on)


% check arguments
if (nargin < 5)
  error('solveRxnDiffusionEqn1dForwardEuler: missing arguments');
end
if (nargin < 6)
  debug_on = 0;
end
if (nargin < 7)
  timing_on = 0;
end

% construct grid
x_lo = -10;
x_hi = 10;
dx = (x_hi-x_lo)/N;
x = x_lo:dx:x_hi;  x = x';

% construct Laplacian operator (with boundary conditions) 
e = ones(N+1,1);                         
L = 1/dx^2*spdiags([e -2*e e], -1:1, N+1, N+1);
L(1,1) = 0; L(1,2) = 0;       % no need to compute Laplacian at x = -10
L(N+1,N) = 0; L(N+1,N+1) = 0; % no need to compute Laplacian at x =  10

% start clock to measure computation 
if (timing_on == 1)
  t_start = cputime;
end

% set initial conditions
u = 1./(1+x.^2);
v = x./(1+x.^2);

% forward Euler time integration
t = t_init;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = 1./(1+(x-t).^2);
    v_exact = (x-t)./(1+(x-t).^2);
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
    set(gca, 'XLim', [x_lo x_hi]);
    title_string = sprintf('t = %f',t);
    title(title_string);

    % plot current error
    figure(2); clf;
    plot(x,err_u,'b');
    hold on;
    plot(x,err_v,'g');
    set(gca, 'XLim', [x_lo x_hi]);
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow
    pause

  end %  end case: (debug_on == 1)

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t; 
  end

  % impose boundary conditions
  u(1)   = 1/(1+(x_lo-t)^2);
  u(end) = 1/(1+(x_hi-t)^2);
  v(1)   = (x_lo-t)/(1+(x_lo-t)^2);
  v(end) = (x_hi-t)/(1+(x_hi-t)^2);

  % compute u_t and v_t 
  u_t = L*u + 2*(v.^3 + v.*u.^2 - 3*v.^2.*u + u.^3);
  v_t = D*L*v - u.^3 + (x-t).*v.^3 + 6*D*u.^2.*v - 2*D*v.^3;

  % update solution
  u = u + dt*u_t;
  v = v + dt*v_t;

  % update time
  t = t + dt;

  % fix boundary conditions
  u(1)   = 1/(1+(x_lo-t)^2);
  u(end) = 1/(1+(x_hi-t)^2);
  v(1)   = (x_lo-t)/(1+(x_lo-t)^2);
  v(end) = (x_hi-t)/(1+(x_hi-t)^2);

end

% measure time to solve reaction-diffusion equations
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
u_exact = 1./(1+(x-t).^2);
v_exact = (x-t)./(1+(x-t).^2);
