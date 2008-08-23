%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveWaveEqn1dKPY_Transformed_OTS() computes the solutions of the 1d wave 
% equation 
%
%   u_tt = (c(x))^2 u_xx + f(x,t)
%
% on the domain -1 < x < 1 subject to the initial conditions
%
%   u(x,0)   = sin(2*pi*x) + cos(3*pi*x)
%   u_t(x,0) = -pi*(2*cos(2*pi*x) + 3*sin(3*pi*x))
%
% and boundary conditions
%
%   u(-1,t) = sin(2*pi*(-1-t)) + cos(3*pi*(-1+t))
%   u(1,t)  = sin(2*pi*(1-t)) + cos(3*pi*(1+t))
%
% The wave speed and source term are given by
%
%   c(x) = 1 + 0.5*sin(0.5*pi*x)
%
%   f(x,t) = 0.25*pi^2*sin(0.5*pi*x) ...
%          * ( 16*sin(2*pi*(x-t)) + 36*cos(3*pi*(x+t)) ...
%            + 4*sin(0.5*pi*x)*sin(2*pi*(x-t)) ...
%            + 9*sin(0.5*pi*x)*cos(3*pi*(x+t)) )
%
% The analytical solution to this problem is 
%
%   u(x,t) = sin(2*pi*(x-t)) + cos(3*pi*(x+t))
%
% We solve the variable coefficient wave equation on the transformed domain 
% 
%   y = 4/pi*( atan((2*tan(pi*x/4)+1)/sqrt(3)) + pi/6 ) - 1
%
% so that the leading-order spatial derivative has a constant coefficient.
% The wave equation in the transformed domain is given by
%
%   u_tt = 4 c_bar^2 u_yy - 2 c_bar c'(x) u_y  + f(x,t)
%
% on the domain -1 < y < 1 with boundary conditions 
%
%   u(-1,t) = sin(2*pi*(-1-t)) + cos(3*pi*(-1+t))
%   u(1,t)  = sin(2*pi*(1-t)) + cos(3*pi*(1+t))
%
% appropriately transformed initial conditions.
%
% The numerical solution is computed a node-centered grid using the
% direct completely centered discretization of the second-order wave 
% equation introduced by Kreiss, Petersson, and Ystrom (2002).  The optimal
% time step of dt = dy/(2*c_bar) is used and correction terms for the 
% lower-order spatial derivative and source terms are used.
%
% USAGE:
%   function [u, u_exact, x] = solveWaveEqn1dKPY_Transformed_OTS(N, ...
%                                                                t_final, ...
%                                                                debug_on)
%
% Arguments:
% - N:                   number of grid cells to use for computation.
%                        NOTE:  num grid points = (num grid cells + 1)
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
% 2008/08:  Initial version of code. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x] = solveWaveEqn1dKPY_Transformed_OTS(N, ...
                                                             t_final, ...
                                                             debug_on)


% check arguments
if (nargin < 2)
  error('solveWaveEqn1dKPY_Transformed_OTS: missing arguments');
end
if (nargin < 3)
  debug_on = 0;
end

% construct grid
y_lo = -1.0;
y_hi =  1.0;
dy = (y_hi-y_lo)/N;
y = y_lo:dy:y_hi; y = y';  

e = ones(N+1,1);
L = 1/dy^2*spdiags([e -2*e e], -1:1, N+1, N+1);
L(1,:) = 0;    % no need to update boundary condition
L(end,:) = 0;  % no need to update boundary condition

% construct gradient operator (with boundary conditions) 
e = ones(N+1,1);
G = 0.5/dy*spdiags([e -e], [1,-1], N+1, N+1);
G(1,:) = 0;    % no need to update boundary condition
G(end,:) = 0;  % no need to update boundary condition

% compute transformation from x to y
x = 4/pi*atan(0.5*(sqrt(3)*tan(pi/4*(y+1)-pi/6) - 1));

% set c_bar
c_bar = sqrt(3)/4;

% set c_prime
c_prime = 0.25*pi*cos(0.5*pi*x);

% compute optimal time step
dt = 0.5*dy/c_bar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KPY time integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize t
t = 0.0;

% set initial conditions
u   = sin(2*pi*x) + cos(3*pi*x);
u_t = -pi*(2*cos(2*pi*x) + 3*sin(3*pi*x));

% use second-order forward Euler for first time step
f =  0.25*pi^2*sin(0.5*pi*x) ...
  .* ( 16*sin(2*pi*x) + 36*cos(3*pi*x) ...
     + 4*sin(0.5*pi*x).*sin(2*pi*x)   ...
     + 9*sin(0.5*pi*x).*cos(3*pi*x) );
f_t = -0.25*pi^3*sin(0.5*pi*x) ...
    .* ( 32*cos(2*pi*x) + 108*sin(3*pi*x) ...
       + 8*sin(0.5*pi*x).*cos(2*pi*x) ...
       + 27*sin(0.5*pi*x).*sin(3*pi*x) );
u_tt = 4*c_bar^2*(L*u) - 2*c_bar*c_prime.*(G*u) + f;
u_ttt = 4*c_bar^2*(L*u_t) -2*c_bar*c_prime.*(G*u_t) + f_t;
u_next = u + dt*u_t + 0.5*dt^2*u_tt + 1/6*dt^3*u_ttt;

% update u_prev and u
u_prev = u;
u = u_next;

% update time
t = t + dt;

% update boundary conditions
u(1) = sin(2*pi*(-1-t)) + cos(3*pi*(-1+t));
u(end) = sin(2*pi*(1-t)) + cos(3*pi*(1+t));

while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = sin(2*pi*(x-t)) + cos(3*pi*(x+t));
    err = u-u_exact;
    err_L_inf = norm(err,'inf')

    % plot current solution
    figure(1); clf;
    plot(x,u,'bo');
    hold on;
    plot(x,u_exact,'r');
    title_string = sprintf('t = %f',t);
    axis([y_lo y_hi -5 5]);
    title(title_string);

    % plot current error
    figure(2); clf;
    plot(x,err);
    title_string = sprintf('t = %f',t);
    title(title_string);
    drawnow
    pause

  end %  end case: (debug_on == 1)

  % compute source term
  f =  0.25*pi^2*sin(0.5*pi*x) ...
    .* ( 16*sin(2*pi*(x-t)) + 36*cos(3*pi*(x+t)) ...
       + 4*sin(0.5*pi*x).*sin(2*pi*(x-t)) ...
       + 9*sin(0.5*pi*x).*cos(3*pi*(x+t)) );

  % compute u_tt
  u_tt = 4*c_bar^2*(L*u) - 2*c_bar*c_prime.*(G*u) + f;
  
  % compute u_tt at end points from boundary condition
  u_tt(1)   = -4*pi^2*sin(2*pi*(-1-t)) - 9*pi^2*cos(3*pi*(-1+t));
  u_tt(end) = -4*pi^2*sin(2*pi*(1-t)) - 9*pi^2*cos(3*pi*(1+t));

  % compute correction term for u_tt
  f_tt = -0.25*pi^4*sin(0.5*pi*x) ...
       .*( 64*sin(2*pi*(x-t)) + 324*cos(3*pi*(x+t)) ...
         + 16*sin(0.5*pi*x).*sin(2*pi*(x-t)) ...
         + 81*sin(0.5*pi*x).*cos(3*pi*(x+t)) );
  u_tt_corr = dt^2/12*(-16*c_bar^3*(G*c_prime).*(L*u) ...
                      - 8*c_bar^3*(L*c_prime).*(G*u) ...
                      + 4*c_bar^2*(c_prime.^2).*(L*u) ...
                      + 4*c_bar^2*c_prime.*(G*c_prime).*(G*u) ...
                      + 4*c_bar^2*L*f - 2*c_bar*c_prime.*(G*f) ...
                      + f_tt);

  % update solution
  u_next = 2*u - u_prev + dt^2*(u_tt + u_tt_corr);

  % update time
  t = t + dt;

  % update u_prev and u 
  if (t < t_final)
    u_prev = u;
    u = u_next;

  else
    % if we have overstepped t_final, use quadratic interpolation to obtain
    % u(t_final).
    % NOTE:  solution is third-order accurate because error in quadratic 
    %        interpolation is O(dt^3).
    dt_final = t_final - (t-dt);
    u = 1/dt^2*( 0.5*dt_final*(dt_final-dt)*u_prev ...
               - (dt_final-dt)*(dt_final+dt)*u ...
               + 0.5*dt_final*(dt_final+dt)*u_next );
    t = t_final;
  end

  % update boundary conditions
  u(1) = sin(2*pi*(-1-t)) + cos(3*pi*(-1+t));
  u(end) = sin(2*pi*(1-t)) + cos(3*pi*(1+t));

end

% compute exact solution 
u_exact = sin(2*pi*(x-t)) + cos(3*pi*(x+t));
