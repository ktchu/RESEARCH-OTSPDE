%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test_1st_derivative_stencil_5th_order.m checks that the fifth-order 
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
N = [25 50 100 200 400];
stencil = [-1/30, 1/4, -1, 1/3, 1/2, -1/20];  % (-3:2)


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
  for j = 4:length(x)-2
    d1f_approx(j) = sum(f(j-3:j+2).*stencil)/dx;
  end
  err_d1f(i) = norm(d1f_approx(4:end-2)-d1f_exact(4:end-2),'inf');

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
N_plot = 10:10:1000;
plot(N_plot,exp(log(N_plot)*P(1)+P(2)),'r');
order_str = sprintf('Order = %1.1f', order);
text(100,0.001,order_str);
