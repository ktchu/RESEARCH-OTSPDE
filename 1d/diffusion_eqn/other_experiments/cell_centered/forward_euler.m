%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script computes the solutions of the 1d heat equation
%
%   T_t = kappa T_xx + f
%
% on the domain 0 < x < 1 with a Dirichlet boundary condition imposed at
% x = 0 and a Neumann boundary condition imposed at x = 1.  The numerical
% solution is computed on a cell-centered grid using forward Euler time
% integration with a second-order central difference approximation to the
% Laplacian.  A suboptimal time step is used.  
%
% The discretization for the boundary condition can be set to 
% linear extrapolation or quadratic extrapolation.
%  
% Kevin T. Chu
% 2007 September
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% physical parameters
kappa = 2.0;  % thermal conductivity

% boundary conditions parameters
% NOTE: BC_order should be set to 1 for linear extrapolation and
%       2 for quadratic extrapolation.  An error will result for
%       any other values.
BC_order = 1;
T_0 = 1;
dTdx_1 = 0.5;

% time integration parameters
t_init  = 0.0;
t_final = 0.01;

% construct grid
N = 100;
dx = 1/N;
x = -dx/2:dx:1+dx/2;  x = x';   % include ghost cells for boundary conditions

% set dt
dt = dx^2/kappa/3;

% construct Laplacian operator (with boundary conditions) 
L = (diag(ones(N+1,1),1) - 2*diag(ones(N+2,1),0) + diag(ones(N+1,1),-1))/dx/dx;
L(1,1) = 0; L(1,2) = 0;   % no need to compute Laplacian at Dirichlet boundary 
L(N+2,N+1) = 0; L(N+2,N+2) = 0;  % no need to compute Laplacian for 
                                 % ghost cell at Neumann boundary

% initialize temperature
T = 1.0 + 0.5*x + 2*sin(5/2*pi*x) - 4*sin(11/2*pi*x) + 3*sin(7/2*pi*x);
f = 5*sin(3/2*pi*x) - 7*sin(15/2*pi*x) + 10*sin(21/2*pi*x);
f = 0;


% forward Euler time integration
t = t_init;
figure(1); clf; % clear figure
while (t < t_final)
  
  % compute exact solution
  if (f == 0)
    T_exact = 1.0 + 0.5*x ...
            + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*kappa*t) ...
            - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*kappa*t) ...
            + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*kappa*t);
  else 
    T_exact = 1.0 + 0.5*x ...
            + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*kappa*t) ...
            - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*kappa*t) ...
            + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*kappa*t) ...
            + 5/kappa*4/9/pi^2*sin(3/2*pi*x)*(1-exp(-9/4*pi^2*kappa*t)) ...
            - 7/kappa*4/225/pi^2*sin(15/2*pi*x)*(1-exp(-225/4*pi^2*kappa*t)) ...
            + 10/kappa*4/21^2/pi^2*sin(21/2*pi*x)*(1-exp(-21^2/4*pi^2*kappa*t));
  end
  err = T(2:N+1)-T_exact(2:N+1);
  err_L_inf = norm(err,'inf')

  % plot current solution
  figure(1); clf;
  plot(x,T,'bo');
  hold on;
  plot(x,T_exact,'r');
  plot_axes = axis; plot_axes(1) = 0; plot_axes(2) = 1;
  axis(plot_axes);
  title_string = sprintf('t = %f',t);
  title(title_string);

  % plot current error
  figure(2); clf;
  plot(x(2:N+1),err,'b');
  plot_axes = axis; plot_axes(1) = 0; plot_axes(2) = 1;
  axis(plot_axes);
  title_string = sprintf('t = %f',t);
  title(title_string);
  drawnow

  % adjust dt so we don't overstep t_final
  if (t+dt > t_final)
    dt = t_final - t;
  end

  % impose boundary conditions
  if (BC_order == 1)
    T(1) = 2*T_0 - T(2);
  elseif (BC_order == 2)
    T(1) = 8/3*T_0 - 2*T(2) + T(3)/3;
  else
    error('Invalid BC_order');
  end
  T(N+2) = T(N+1) + dx*dTdx_1;

  % update solution
  if (f == 0)
    T = T + kappa*dt*(L*T);
  else 
    T = T + kappa*dt*(L*T) + dt*f; 
  end

  % update time
  t = t + dt;

end

% compute exact solution and error
if (f == 0)
  T_exact = 1.0 + 0.5*x ...
          + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*kappa*t) ...
          - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*kappa*t) ...
          + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*kappa*t);
else
  T_exact = 1.0 + 0.5*x ...
          + 2*sin(5/2*pi*x)*exp(-25/4*pi^2*kappa*t) ...
          - 4*sin(11/2*pi*x)*exp(-121/4*pi^2*kappa*t) ...
          + 3*sin(7/2*pi*x)*exp(-49/4*pi^2*kappa*t) ...
          + 5/kappa*4/9/pi^2*sin(3/2*pi*x)*(1-exp(-9/4*pi^2*kappa*t)) ...
          - 7/kappa*4/225/pi^2*sin(15/2*pi*x)*(1-exp(-225/4*pi^2*kappa*t)) ...
          + 10/kappa*4/21^2/pi^2*sin(21/2*pi*x)*(1-exp(-21^2/4*pi^2*kappa*t));
end
err = T(2:N+1)-T_exact(2:N+1);
err_L_inf = norm(err,'inf')


% impose boundary conditions
if (BC_order == 1)
  T(1) = 2*T_0 - T(2);
elseif (BC_order == 2)
  T(1) = 8/3*T_0 - 2*T(2) + T(3)/3;
else
  error('Invalid BC_order');
end
T(N+2) = T(N+1) + dx*dTdx_1;

% plot final result
figure(1); clf;
plot(x,T,'bo')
hold on;
plot(x,T_exact,'r')
plot_axes = axis; plot_axes(1) = 0; plot_axes(2) = 1;
axis(plot_axes);
title_string = sprintf('t = %f',t);
title(title_string);

% plot error
figure(2); clf;
plot(x(2:N+1),err,'b');
plot_axes = axis; plot_axes(1) = 0; plot_axes(2) = 1;
axis(plot_axes);
title_string = sprintf('t = %f',t);
title(title_string);
err2 = abs(T(2)-T_exact(2))
