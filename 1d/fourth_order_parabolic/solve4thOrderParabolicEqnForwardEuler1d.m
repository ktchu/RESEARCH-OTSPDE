%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solve4thOrderParabolicEqnForwardEuler1d() computes the solutions of 
% the 4th-order parabolic equation in one spatial dimensions
%
%   u_t = -u_xxxx + f
%
% where
%
%   f(x,t) = -3*sin(5*pi*x)*sin(3*t)
%          + 625*pi^4*sin(5*pi*x)*cos(3*t)
%          + 4*pi^2*sin(2*pi*x)*exp(-4*pi^2*t)
%          + 16*pi^4*sin(2*pi*x)*(3-exp(-4*pi^2*t))
%
% on the domain -1 < x < 1 with homogeneous Dirichlet boundary conditions 
% imposed at all boundaries.  The numerical solution is computed on a 
% node-centered grid using forward Euler time integration with a 
% fourth-order central difference approximation to the bilaplacian.  
%
% USAGE:
%   function [u, u_exact, x, timing_data] = ...
%     solve4thOrderParabolicEqnForwardEuler1d(dx, dt, ...
%                                             t_init, t_final, ...
%                                             debug_on, timing_on)
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
% - u:                   numerical solution
% - u_exact:             analytical solution
% - x:                   grid points
% - timing_data:         array of the following form containing timing data
%                        [bilaplacian_construction time, solution time].
%                        If timing is not activated, timing_data is
%                        set to [-1, -1].
%
%
% NOTES: 
% - For the purposes of illustration, the boundary conditions for the
%   example are imposed by filling the ghostcells from the exact
%   solution.
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
  solve4thOrderParabolicEqnForwardEuler1d(dx, dt, ...
                                          t_init, t_final, ...
                                          debug_on, timing_on)


% check arguments
if (nargin < 4)
  error('solve4thOrderParabolicEqnForwardEuler1d: missing arguments');
end
if (nargin < 5)
  debug_on = 0;
end
if (nargin < 6)
  timing_on = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for main computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  
% construct grid
N = 2/dx + 5;
x = -1-2*dx:dx:1+2*dx; x = x'; % include ghostcells

% start clock to measure time to prepare for computation
disp('Constructing Bilaplacian operators ...');
if (timing_on == 1)
  t_start = cputime;
end

% construct 4nd-order Biaplacian operators (with boundary conditions) 
e = ones(N,1);
L = 1/dx^4 * spdiags([-1/6*e 2*e -13/2*e 28/3*e -13/2*e 2*e -1/6*e], ...
                      -3:3, N, N);
L(1,:) = 0; L(2,:) = 0; L(3,:) = 0;
L(N,:) = 0; L(N-1,:) = 0; L(N-2,:) = 0;

% measure Bilaplacian construction time 
if (timing_on == 1)
  time_bilaplacian_construction = cputime - t_start;
end

% set up plotting grid
if (debug_on == 1) 
  % maximum number of grid cells to use for plotting
  N_plot = 101;
  if (N < N_plot)
    N_plot = N;
  end
  dx_plot = 2/(N_plot-1);
  x_plot = -1:dx_plot:1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve fourth-order parabolic equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% restart clock to measure computation time
if (timing_on == 1)
  t_start = cputime;
end

% initialize temperature
u = 2*sin(2*pi*x) + sin(5*pi*x);

% forward Euler time integration
t = t_init;
time_step = 1;
while (t < t_final)
 
  if (debug_on == 1)
 
    % compute exact solution
    u_exact = 2*sin(2*pi*x)*cos(10000*pi*t) + sin(5*pi*x)*(2-exp(-100*pi^2*t));
    err = u-u_exact;
    err(1:3) = 0; err(N-2:N) = 0;
    err_L_inf = norm(err,'inf');

    % plot current solution
    figure(1); clf;
    u_plot = interp1(x, u, x_plot, 'cubic');
    plot(x_plot,u_plot,'bo');
    hold on;
    u_plot = interp1(x, u_exact, x_plot, 'cubic');
    plot(x_plot,u_plot,'r-');
    title_string = sprintf('t = %g',t);
    title(title_string);
    xlabel('x'); 

    % plot current error
    figure(2); clf;
    err_plot = interp1(x, err, x_plot, 'cubic');
    plot(x_plot,err_plot);
    title_string = sprintf('t = %g',t);
    title(title_string);
    xlabel('x'); 
    drawnow

    % display current status of solution
    status_str = sprintf('L infinity Error = %g', err_L_inf);
    disp(status_str); 

  end % end case: (debug_on == 1)

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t;
  end

  % impose boundary conditions
  u(1:3) = 2*sin(2*pi*x(1:3))*cos(10000*pi*t) ...
         + sin(5*pi*x(1:3))*(2-exp(-100*pi^2*t));
  u(N-2:N) = 2*sin(2*pi*x(N-2:N))*cos(10000*pi*t) ...
           + sin(5*pi*x(N-2:N))*(2-exp(-100*pi^2*t));

  % compute source terms
  f = -20000*pi*sin(2*pi*x)*sin(10000*pi*t) ...
    + 32*pi^4*sin(2*pi*x)*cos(10000*pi*t) ...
    + 100*pi^2*sin(5*pi*x)*exp(-100*pi^2*t) ...
    + 625*pi^4*sin(5*pi*x)*(2-exp(-100*pi^2*t));

  % update solution
  u = u + dt*(-L*u + f);

  % update time and time step
  t = t + dt;
  time_step = time_step + 1;

end

% measure time to solve fourth-order parabolic equation
if (timing_on == 1)
  time_solve = cputime - t_start;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output timing statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (timing_on)
  timing_data = [time_bilaplacian_construction, time_solve];

  if (debug_on)
    disp(' ');
    disp('==================================================================');
    disp('Computation Statistics');
    bilapl_construct_time_str = ...
      sprintf('  Bilaplacian Construction Time: %f', ...
              time_bilaplacian_construction);
    comp_time_str = sprintf('  Solution Time: %f', time_solve);
    disp(bilapl_construct_time_str);
    disp(comp_time_str);
    disp('==================================================================');
  end

else
  
  timing_data = [-1, -1];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute exact solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_exact = 2*sin(2*pi*x)*cos(10000*pi*t) + sin(5*pi*x)*(2-exp(-100*pi^2*t));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove ghost cells from numerical soln
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = u(3:N-2);
u_exact = u_exact(3:N-2);
x = x(3:N-2);
