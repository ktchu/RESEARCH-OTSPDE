%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveRxnDiffusionEqn1dForwardEuler() computes the solutions of the 
% following 1d reaction-diffusion equation 
%
%   u_t = u_xx - u*(u-2)*(u+1)
%
% on the domain -5 < x < 5 with boundary conditions imposed from the 
% exact travelling-wave solution:  
%
%   u(x,t) = 2*(exp(eta) - 1)/(exp(eta) + 2),
%
% where eta = (k*x + w*t) with k = 3/sqrt(2) and w = 3/2, which implies a 
% wave speed of -1/sqrt(2).  The initial conditions at t = 0 is given by 
%
%   u(x,0) = 2*(exp(k*x) - 1)/(exp(k*x) + 2).
%
% The numerical solution is computed on a node-centered grid using 
% forward Euler time integration with a second-order central difference 
% approximation to the Laplacian.  
%
% USAGE:
%   function [u, u_exact, x, timing_data] = ...
%     solveRxnDiffusionEqn1dForwardEuler(dx, dt, ...
%                                        t_final, ...
%                                        debug_on, timing_on)
%
% Arguments:
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
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x, timing_data] = ...
  solveRxnDiffusionEqn1dForwardEuler(dx, dt, ...
                                     t_final, ...
                                     debug_on, timing_on)


% check arguments
if (nargin < 3)
  error('solveRxnDiffusionEqn1dForwardEuler: missing arguments');
end
if (nargin < 4)
  debug_on = 0;
end
if (nargin < 5)
  timing_on = 0;
end

% construct grid
x_lo = -5;
x_hi =  5;
x = x_lo:dx:x_hi;  x = x';   
N = length(x);

% compute k and w
k = 3/sqrt(2);
w = 3/2;

% construct Laplacian operator (with boundary conditions) 
e = ones(N,1);                         
L = 1/dx^2*spdiags([e -2*e e], -1:1, N, N);
L(1,1) = 0; L(1,2) = 0;  % no need to compute Laplacian at boundary 
L(N,N) = 0; L(N,N) = 0;  % no need to compute Laplacian at boundary

% start clock to measure computation 
if (timing_on == 1)
  t_start = cputime;
end

% set initial conditions
u = 2*(exp(k*x) - 1)./(exp(k*x) + 2);

% forward Euler time integration
t = 0.0;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = 2*(exp(k*x+w*t) - 1)./(exp(k*x+w*t) + 2);
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
  u = u + dt*(L*u - u.*(u-2).*(u+1));

  % update time
  t = t + dt;

  % impose boundary conditions
  u(1) = 2*(exp(k*x_lo+w*t) - 1)/(exp(k*x_lo+w*t) + 2);
  u(end) = 2*(exp(k*x_hi+w*t) - 1)/(exp(k*x_hi+w*t) + 2);

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
    comp_t_solve = sprintf('  Solution Time: %f', t_solve);
    disp(comp_time_str);
    disp('==================================================================');
  end

else

  timing_data = -1;

end

% compute exact solution 
u_exact = 2*(exp(k*x+w*t) - 1)./(exp(k*x+w*t) + 2);
