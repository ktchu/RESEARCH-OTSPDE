%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveDiffusionEqn1dDuFortFrankelOTS() computes the solutions of the 1d 
% diffusion equation 
%
%   u_t = D u_xx + f(x)
%
% on the domain 0 < x < 1 with a Dirichlet boundary condition imposed at
% x = 0 and a Neumann boundary condition imposed at x = 1.  The numerical
% solution is computed on a node-centered grid using the DuFort-Frankel
% scheme.  The optimal time step of dx^2/(sqrt(12) D) and correction terms
% are used.
%
% USAGE:
%   function [u, u_exact, x, timing_data] = ...
%     solveDiffusionEqn1dDuFortFrankelOTS(D, ...
%                                         source_term_type, ...
%                                         u_0, dudx_1, ...
%                                         dx, ...
%                                         t_init, t_final, ...
%                                         debug_on, timing_on)
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
% 2009/05:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x, timing_data] = ...
  solveDiffusionEqn1dDuFortFrankelOTS(D, ...
                                      source_term_type, ...
                                      u_0, dudx_1, ...
                                      dx, ...
                                      t_init, t_final, ...
                                      debug_on, timing_on)


% check arguments
if (nargin < 7)
  error('solveDiffusionEqn1dDuFortFrankelOTS: missing arguments');
end
if (nargin < 8)
  debug_on = 0;
end
if (nargin < 9)
  timing_on = 0;
end

% construct grid
N = 1/dx + 1;
x = 0:dx:1+dx;  x = x';   % include ghost cell for Neumann BC

% compute optimal time step
dt = dx^2/sqrt(12)/D;

% construct Laplacian operator (with boundary conditions) 
e = ones(N+1,1);                         
L = 1/dx^2*spdiags([e -2*e e], -1:1, N+1, N+1);
L(1,1) = 0; L(1,2) = 0;   % no need to compute Laplacian at Dirichlet boundary 
L(N+1,N) = 0; L(N+1,N+1) = 0;  % no need to compute Laplacian for 
                               % ghost cell at Neumann boundary

% cache sigma = dt/dx^2
sigma = dt/dx^2;

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
  error('solveDiffusionEqn1dDuFortFrankelOTS: Invalid source term type.  Valid valuesare 0 and 1');
end

% take first time step using forward Euler time integration with OTS
t = t_init;

% impose boundary conditions
u(1) = u_0;
u(N+1) = u(N-1) + 2*dx*dudx_1;

% compute solution at next time step
u_next = u + dt*D*(L*u) + dt*f;

% update time
t = t + dt;

% update solutions
u_prev = u;
u = u_next;


% DuFort-Frankel scheme
while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
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
    plot(x(1:N),err);
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow

  end %  end case: (debug_on == 1)

  % impose boundary conditions
  u(1) = u_0;
  u(N+1) = u(N-1) + 2*dx*dudx_1;

  % compute solution at next time step
  corr_term = 2*D^2*dt^2*sigma*L*f;
  u_next = 1/(1+2*D*sigma) * ( (1-2*D*sigma)*u_prev + 4*D*sigma*u ...
                             + 2*dt*D*(L*u) + 2*dt*f + corr_term ); 

  % update time
  t = t + dt;

  % update solutions
  u_prev = u;
  u = u_next;

  % if we have overstepped t_final, use linear interpolation to get
  % the solution at the desired time
  if (t > t_final)
    alpha = (t - t_final)/dt; 
    u = (1-alpha)*u + alpha*u_prev;
    t = t_final;
  end

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

% remove ghost cell
u = u(1:N);
u_exact = u_exact(1:N);
x = x(1:N);
