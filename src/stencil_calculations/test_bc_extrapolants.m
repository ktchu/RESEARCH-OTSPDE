%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test_bc_extrapolants.m checks that the extrapolants used to set the 
% ghostcells for a Neumann boundary condition in 1d have the expected
% order of accuracy.  Only cubic and quartic extrapolants can be tested.
%
% Set the polynomial_order variable to choose the polynomial extrapolant
% to test.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
format long;
format compact;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
polynomial_order = 4;  % choices: 3, 4
grid_sizes = [100 200 400 800 1600];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute error for different dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_f_extrap = zeros(size(grid_sizes));
for i = 1:length(grid_sizes)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute grid 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  N = grid_sizes(i);
  dx = 1/N;
  x = 0:dx:1+dx;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compare results with analytical derivative
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  f = exp(cos(3*x))+sin(5*x)+x.^8;
  dfdx = -3*sin(3*x).*exp(cos(3*x)) + 5*cos(5*x) + 8*x.^7;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute extrapolated value at x_{N+1}
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if polynomial_order == 3
    % cubic extrapolant
    alpha = 0.5*(dx*dfdx(N) - 1.5*f(N) + 2*f(N-1) - 0.5*f(N-2));
    f_extrap = 3*f(N) - 3*f(N-1) + f(N-2) + 6*alpha; 

  else

    % quartic extrapolant
    alpha = 1/6*(dx*dfdx(N) - 11/6*f(N) + 3*f(N-1) - 1.5*f(N-2) + 1/3*f(N-3));
    f_extrap = 4*f(N) - 6*f(N-1) + 4*f(N-2) - f(N-3) + 24*alpha;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute error
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  err_f_extrap(i) = abs(f_extrap - f(N+1));

end

% show error
err_f_extrap = err_f_extrap 

% compute order of discretization
P = polyfit(log(grid_sizes),log(err_f_extrap),1);
order = -P(1)

% plot error vs grid size
figure(2); clf;
loglog(grid_sizes,err_f_extrap,'bo');
hold on;
N_plot = 100:10:10000;
plot(N_plot,exp(log(N_plot)*P(1)+P(2)),'r');
order_str = sprintf('Order = %1.1f', order);
if polynomial_order == 3
  text(500,1e-12,order_str);
  axis([100 10000 1e-14 1e-4]);
else
  text(300,1e-13,order_str);
  axis([100 10000 1e-15 1e-6]);
end
