%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqnForwardEulerOTS1d() computes the solutions of the 1d 
% diffusion equation 
%
%   u_t = D u_xx + f(x)
%
% on the domain 0 < x < 1 with f(x) = 2*t - 42*x^5 + 6*x and boundary
% conditions given by u(0) = 0 and u'(1) = 4.  The initial condition
% is taken to be u(x,0) = x^7 - x^3.  The numerical solution is computed on 
% a node-centered grid using forward Euler time integration with a 
% second-order central difference approximation to the Laplacian.  The 
% optimal time step of dx^2/(6D) and correction terms are used.
%
% USAGE:
%   function [u, u_exact, x, timing_data] = ...
%     solveDiffusionEqnForwardEulerOTS1d(D, ...
%                                        dx, ...
%                                        t_final, ...
%                                        debug_on, timing_on)
%
% Arguments:
% - D:                   diffusion coefficient 
% - dx:                  grid spacing
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
  solveDiffusionEqnForwardEulerOTS1d(D, ...
                                     dx, ...
                                     t_final, ...
                                     debug_on, timing_on)


% check arguments
if (nargin < 3)
  error('solveDiffusionEqnForwardEulerOTS1d: missing arguments');
end
if (nargin < 4)
  debug_on = 0;
end
if (nargin < 5)
  timing_on = 0;
end

% construct grid
N = 1/dx + 1;
x = 0:dx:1+dx;  x = x';   % include ghost cell for Neumann BC

% compute optimal time step
dt = dx^2/6/D;

% construct Laplacian operator (with boundary conditions) 
e = ones(N+1,1);                         
L = 1/dx^2*spdiags([e -2*e e], -1:1, N+1, N+1);
L(1,1) = 0; L(1,2) = 0;   % no need to compute Laplacian at Dirichlet boundary 
L(N+1,N) = 0; L(N+1,N+1) = 0;  % no need to compute Laplacian for 
                               % ghost cell at Neumann boundary

% start clock to measure computation 
if (timing_on == 1)
  t_start = cputime;
end

% set initial conditions
u = sinh(sin(4*pi*x)+cos(3*pi*x)-1);

% forward Euler time integration
t = 0;
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = cos(5*t)*sinh(sin(4*pi*x)+cos(3*pi*x)-1);
    err = u(1:N)-u_exact(1:N);
    err_L_inf = norm(err,'inf')

    % plot current solution
    figure(1); clf;
    plot(x(1:N),u(1:N),'bo');
    hold on;
    plot(x(1:N),u_exact(1:N),'r');
    title_string = sprintf('t = %f',t);
    title(title_string);

    % plot current error
    figure(2); clf;
    plot(x(1:N),err,'o');
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow

    pause

  end %  end case: (debug_on == 1)

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t; 
  end

  % set boundary conditions
  u_0 = 1e6*dx^4;
  u_0 = 0;
  dudx_1 = 4*pi*cos(5*t)*cosh(2);

  % set source terms
  f = -5*sin(5*t)*sinh(sin(4*pi*x)+cos(3*pi*x)-1) ...
    - D*cos(5*t)*(4*pi*cos(4*pi*x)-3*pi*sin(3*pi*x)).^2 ...
                .*sinh(sin(4*pi*x)+cos(3*pi*x)-1) ...
    - D*cos(5*t)*(-16*pi^2*sin(4*pi*x)-9*pi^2*cos(3*pi*x)) ...
                .*cosh(sin(4*pi*x)+cos(3*pi*x)-1);
  f_t = -25*cos(5*t)*sinh(sin(4*pi*x)+cos(3*pi*x)-1) ...
      + 5*D*sin(5*t)*((4*pi*cos(4*pi*x)-3*pi*sin(3*pi*x)).^2 ...
                     .*sinh(sin(4*pi*x)+cos(3*pi*x)-1) ...
                     +(-16*pi^2*sin(4*pi*x)-9*pi^2*cos(3*pi*x)) ...
                     .*cosh(sin(4*pi*x)+cos(3*pi*x)-1));

  % impose boundary conditions
  u(1) = u_0;

  % 3rd-order ghostcell value
  %u(N+1) = u(N-1) + 2*dx*dudx_1;  
  %u(N+1) = 2*u(N) - 2.5*u(N-1) + 2*u(N-2) - 0.5*u(N-3) + dx*dudx_1;

  % 4th-order ghostcell value
  %u(N+1) = -1.5*u(N) + 3*u(N-1) - 0.5*u(N-2) + 3*dx*dudx_1;  
  %u(N+1) = 1.5*u(N) - u(N-1) + 0.5*u(N-2) + dx*dudx_1;  
  %u(N+1) = 13/6*u(N) - 3*u(N-1) + 2.5*u(N-2) - 2/3*u(N-3) + dx*dudx_1;
  %alpha = 0.5*(dx*dudx_1 - 1.5*u(N) + 2*u(N-1) - 0.5*u(N-2));
  %u(N+1) = 3*u(N) - 3*u(N-1) + u(N-2) + 6*alpha;

  % 5th-order ghostcell value
  %u(N+1) = -10/3*u(N) + 6*u(N-1) - 2*u(N-2) + 1/3*u(N-3) + 4*dx*dudx_1;
  %u(N+1) = 5/6*u(N) - 2*u(N-1) + 4*u(N-2) - 7/3*u(N-3) + 0.5*u(N-4) ...
  %       + 2*dx*dudx_1;
  %u(N+1) = -1.25*u(N) + 2*u(N-1) + u(N-2) - u(N-3) + 0.25*u(N-4) ...
  %       + 3*dx*dudx_1;
  %u(N+1) = 13/6*u(N) - 6*u(N-1) + 7*u(N-2) - 11/3*u(N-3) + 0.75*u(N-4) ...
  %       + dx*dudx_1;
  %alpha = 1/6*(dx*dudx_1 - 11/6*u(N) + 3*u(N-1) - 1.5*u(N-2) + 1/3*u(N-3)); 
  %u(N+1) = 4*u(N) - 6*u(N-1) + 4*u(N-2) - u(N-3) + 24*alpha; % - 2e5*dx^5;

  % 6th-order ghostcell value
  u(N+1) = -65/12*u(N) + 10*u(N-1) - 5*u(N-2) + 5/3*u(N-3) - 0.25*u(N-4) ...
         + 5*dx*dudx_1;

  % exact ghostcell value
  %u(N+1) - cos(5*t)*sinh(sin(4*pi*x(end))+cos(3*pi*x(end))-1)
  %u(N+1) = cos(5*t)*sinh(sin(4*pi*x(end))+cos(3*pi*x(end))-1);

  % update solution
  u = u + dt*D*(L*u) + dt*f + 0.5*dt^2*(D*L*f + f_t);
  u(1) = u_0;

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

% remove ghost cell
u = u(1:N);
u_exact = u_exact(1:N);
x = x(1:N);
