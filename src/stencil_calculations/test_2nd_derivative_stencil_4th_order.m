%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test_2nd_derivative_stencil_4th_order.m checks that the fourth-order 
% accurate stencil for the second derivative operator is correct and 
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
N = [50 100 200 400 800 1000];
stencil = [-1/12, 4/3, -5/2, 4/3, -1/12];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute error for different dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_d2f = zeros(size(N));
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
  d2f_exact = -9*cos(3*x).*exp(cos(3*x)) ...
            + 9*(sin(3*x).^2).*exp(cos(3*x)) ...
            - 25*sin(5*x) + 56*x.^6;
  d2f_approx = zeros(size(x));
  for j = 3:length(x)-2
    d2f_approx(j) = sum(f(j-2:j+2).*stencil)/dx^2;
  end
  err_d2f(i) = norm(d2f_approx(4:end-3)-d2f_exact(4:end-3),'inf');

  % plot results
  figure(1); clf;
  plot(x,d2f_exact,'b');
  hold on;
  plot(x(4:end-3),d2f_approx(4:end-3),'rx');

end

% compute order of discretization
P = polyfit(log(N),log(err_d2f),1);
order = -P(1)

% plot error vs N
figure(2); clf;
loglog(N,err_d2f,'bo');
hold on;
N_plot = 10:10:10000;
plot(N_plot,exp(log(N_plot)*P(1)+P(2)),'r');
order_str = sprintf('Order = %1.1f', order);
text(100,1,order_str);
