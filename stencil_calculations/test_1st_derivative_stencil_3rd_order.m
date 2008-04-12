%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test_1st_derivative_stencil_3rd_order.m checks that the third-order 
% accurate stencil for the first derivative operator is correct and 
% yields the expected order of accuracy.
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
N = [100 200 400 800 1000 2000];
stencil = [-1/3, -1/2, 1, -1/6];  % partial upwind (-1:2)
stencil = [1/6, -1, 1/2, 1/3];  % partial upwind (-2:1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute error for different dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_d1f = zeros(size(N));
for i = 1:length(N)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compute grid and stencil
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dx = 2/N(i);
  x = -1:dx:1;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % compare results of applying stencil to known function 
  % with analytical formulae
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  f = exp(cos(3*x))+sin(5*x)+x.^8;
  d1f_exact = -3*sin(3*x).*exp(cos(3*x)) + 5*cos(5*x) + 8*x.^7;
  d1f_approx = zeros(size(x));

  % for partial upwind (-1:2)
  for j = 2:length(x)-2
    d1f_approx(j) = sum(f(j-1:j+2).*stencil)/dx;
  end
  err_d1f(i) = norm(d1f_approx(2:end-2)-d1f_exact(2:end-2),'inf');

  % for partial upwind (-2:1)
  for j = 3:length(x)-1
    d1f_approx(j) = sum(f(j-2:j+1).*stencil)/dx;
  end
  err_d1f(i) = norm(d1f_approx(4:end-3)-d1f_exact(4:end-3),'inf');

  % plot results
  figure(1); clf;
  plot(x,d1f_exact,'b');
  hold on;
  plot(x(4:end-3),d1f_approx(4:end-3),'rx');

end

% compute order of discretization
P = polyfit(log(N),log(err_d1f),1);
order = -P(1)

% plot error vs N
figure(2); clf;
loglog(N,err_d1f,'bo');
hold on;
N_plot = 100:10:10000;
plot(N_plot,exp(log(N_plot)*P(1)+P(2)),'r');
order_str = sprintf('Order = %1.1f', order);
text(1000,0.001,order_str);
