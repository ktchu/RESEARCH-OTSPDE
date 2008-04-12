%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test_4th_derivative_stencil_2nd_order.m checks that the second-order 
% accurate stencil for the fourth derivative operator is correct and
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
stencil = [1, -4, 6, -4, 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute error for different dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_d4f = zeros(size(N));
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
  d4f_exact = 81*cos(3*x).*exp(cos(3*x)) ...
            - 324*sin(3*x).^2.*exp(cos(3*x)) ...
            + 243*cos(3*x).^2.*exp(cos(3*x)) ...
            - 486*cos(3*x).*sin(3*x).^2.*exp(cos(3*x)) ...
            + 81*sin(3*x).^4.*exp(cos(3*x)) ...
            + 625*sin(5*x) + 1680*x.^4;
  d4f_approx = zeros(size(x));
  for j = 3:length(x)-2
    d4f_approx(j) = sum(f(j-2:j+2).*stencil)/dx^4;
  end
  err_d4f(i) = norm(d4f_approx(4:end-3)-d4f_exact(4:end-3),'inf');

  % plot results
  figure(1); clf;
  plot(x,d4f_exact,'b');
  hold on;
  plot(x(4:end-3),d4f_approx(4:end-3),'rx');

end

% compute order of discretization
P = polyfit(log(N),log(err_d4f),1);
order = -P(1)

% plot error vs N
figure(2); clf;
loglog(N,err_d4f,'bo');
hold on;
N_plot = 10:1:1000;
plot(N_plot,exp(log(N_plot)*P(1)+P(2)),'r');
order_str = sprintf('Order = %1.1f', order);
text(100,100,order_str);
