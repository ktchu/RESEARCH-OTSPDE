%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% solveWaveEqn1dKPY_VarMesh_OTS() computes the solutions of the 1d wave 
% equation 
%
%   u_tt = (c(x))^2 u_xx + f(x,t)
%
% on the domain 0 < x < 1 subject to the initial conditions
%
%   u(x,0)   = sin(2*pi*x) + cos(6*pi*x)
%   u_t(x,0) = -pi*(2*cos(2*pi*x) + 6*sin(6*pi*x))
%
% and boundary conditions
%
%   u(0,t) = sin(2*pi*(-t)) + cos(6*pi*t)
%   u(1,t)  = sin(2*pi*(1-t)) + cos(6*pi*(1+t))
%
% The wave speed and source term are given by
%
%   c(x) = 1 - 0.5*cos(pi*x)
%
%   f(x,t) = pi^2*cos(pi*x) ...
%          * ( -4*sin(2*pi*(x-t)) - 36*cos(6*pi*(x+t)) ...
%            + cos(pi*x)*sin(2*pi*(x-t)) ...
%            + 9*cos(pi*x)*cos(6*pi*(x+t)) )
%
% The analytical solution to this problem is 
%
%   u(x,t) = sin(2*pi*(x-t)) + cos(6*pi*(x+t))
%
% The numerical solution is computed on the variable grid size mesh induced 
% by taking a uniform mesh in a transformed domain where the wave speed is 
% constant.  The time integration scheme is based on the centered 
% discretization of the second-order wave equation introduced by Kreiss, 
% Petersson, and Ystrom (2002).  The optimal time step of dt = dy/c_bar 
% is used.  In this expression, dy is the grid spacing of the uniform mesh
% on the transformed domain and c_bar is the harmonic average of the wave 
% speed.  The expression for c_bar is given below.
%
% On the transformed domain 
%
%   y = 2/pi*( atan(sqrt(3)*tan(pi*x/2)) )
%
% the variable coefficient wave equation above becomes a wave equation
% that has a constant coefficient leading-order spatial derivative term:
%
%   u_tt = c_bar^2 u_yy - c_bar c'(x) u_y  + f(x,t),
%
% where c_bar has a value of sqrt(3)/2.
%
% USAGE:
%   function [u, u_exact, x] = solveWaveEqn1dKPY_VarMesh_OTS(N, ...
%                                                            t_final, ...
%                                                            debug_on)
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
% Kevin T. Chu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, u_exact, x] = solveWaveEqn1dKPY_VarMesh_OTS(N, ...
                                                         t_final, ...
                                                         debug_on)


% check arguments
if (nargin < 2)
  error('solveWaveEqn1dKPY_VarMesh_OTS: missing arguments');
end
if (nargin < 3)
  debug_on = 0;
end

% construct grid
y_lo = 0.0;
y_hi = 1.0;
dy = (y_hi-y_lo)/N;
y = y_lo:dy:y_hi; y = y';  

% compute transformation from y to x
x = 2/pi*atan(tan(0.5*pi*y)/sqrt(3));

% construct Laplacian operator (with boundary conditions) 
diag_plus  = (2./(x(3:end)-x(1:end-2))) .* (1./(x(3:end)-x(2:end-1)));
diag_minus = (2./(x(3:end)-x(1:end-2))) .* (1./(x(2:end-1)-x(1:end-2)));
diag_center = -(diag_plus + diag_minus);
diag_plus = [0; 0; diag_plus];
diag_minus = [diag_minus; 0; 0];
diag_center = [0; diag_center; 0];
L = spdiags([diag_minus diag_center diag_plus], -1:1, N+1, N+1);
L(1,:) = 0;    % no need to update boundary condition
L(end,:) = 0;  % no need to update boundary condition

% compute optimal time step 
c_bar = sqrt(3)/2;
dt = dy/c_bar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KPY time integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize t
t = 0.0;

% set initial conditions
u   = sin(2*pi*x) + cos(6*pi*x);
u_t = -pi*(2*cos(2*pi*x) + 6*sin(6*pi*x));

% compute wave speed
c = 1 - 0.5*cos(pi*x);

% use fifth-order Taylor series expansion for first time step
f = pi^2*cos(pi*x) ...
  .* ( -4*sin(2*pi*(x-t)) - 36*cos(6*pi*(x+t)) ...
     + cos(pi*x).*sin(2*pi*(x-t)) ...
     + 9*cos(pi*x).*cos(6*pi*(x+t)) );
f_t = -2*pi^3*cos(pi*x) ...
    .* (-4*cos(2*pi*(x-t)) - 108*sin(6*pi*(x+t)) ...
       + cos(pi*x).*cos(2*pi*(x-t)) ...
       + 27*cos(pi*x).*sin(6*pi*(x+t)) );
f_tt = -4*pi^4*cos(pi*x) ...
     .*(-4*sin(2*pi*(x-t)) - 324*cos(6*pi*(x+t)) ...
       + cos(pi*x).*sin(2*pi*(x-t)) ...
       + 81*cos(pi*x).*cos(6*pi*(x+t)) );

u_tt = -c.^2.*(4*pi^2*sin(2*pi*x) + 36*pi^2*cos(6*pi*x)) + f;
u_ttt = c.^2.*(L*u_t) + f_t;
u_tttt = c.^2.*(L*u_tt) + f_tt;
u_next = u + dt*u_t + 0.5*dt^2*u_tt + 1/6*dt^3*u_ttt + 1/24*dt^4*u_tttt;

% update u_prev and u
u_prev = u;
u = u_next;

% update time
t = t + dt;

% update boundary conditions
u(1) = sin(2*pi*(-t)) + cos(6*pi*t);
u(end) = sin(2*pi*(1-t)) + cos(6*pi*(1+t));

while (t < t_final)

  if (debug_on == 1)

    % compute exact solution and err
    u_exact = sin(2*pi*(x-t)) + cos(6*pi*(x+t));
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
  f = pi^2*cos(pi*x) ...
    .* ( -4*sin(2*pi*(x-t)) - 36*cos(6*pi*(x+t)) ...
       + cos(pi*x).*sin(2*pi*(x-t)) ...
       + 9*cos(pi*x).*cos(6*pi*(x+t)) );

  % compute u_tt
  u_tt = c.^2.*(L*u) + f;

  % compute correction term for u_tt
  f_tt = -4*pi^4*cos(pi*x) ...
       .*(-4*sin(2*pi*(x-t)) - 324*cos(6*pi*(x+t)) ...
         + cos(pi*x).*sin(2*pi*(x-t)) ...
         + 81*cos(pi*x).*cos(6*pi*(x+t)) );
  u_xx = L*u;
  u_tt_corr = dt^2/12*( (c.^2).*(L*(c.^2).*u_xx) ...
                      + c.^2.*(L*f) + f_tt );

  % update solution
  u_next = 2*u - u_prev + dt^2*(u_tt + u_tt_corr);

  % update time
  t = t + dt;

  % update u_minus_two, u_prev and u 
  if (t < t_final)
    u_minus_two = u_prev;
    u_prev = u;
    u = u_next;

  else
    % if we have overstepped t_final, use cubic interpolation to obtain
    % u(t_final).
    % NOTE:  solution is fourth-order accurate because error in cubic
    %        interpolation is O(dt^4).
    dt_final = t_final - (t-dt);
    u = 1/dt^3*(-1/6*(dt_final+dt)*dt_final*(dt_final-dt)*u_minus_two ...
               + 0.5*(dt_final+2*dt)*dt_final*(dt_final-dt)*u_prev ...
               - 0.5*(dt_final+2*dt)*(dt_final+dt)*(dt_final-dt)*u ...
               + 1/6*(dt_final+2*dt)*(dt_final+dt)*dt_final*u_next );
    t = t_final;
  end

  % update boundary conditions
  u(1) = sin(2*pi*(-t)) + cos(6*pi*t);
  u(end) = sin(2*pi*(1-t)) + cos(6*pi*(1+t));

end

% compute exact solution 
u_exact = sin(2*pi*(x-t)) + cos(6*pi*(x+t));
