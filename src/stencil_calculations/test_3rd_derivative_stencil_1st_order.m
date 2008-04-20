%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test_3rd_derivative_stencil_3rd_order.m checks that the third-order 
% accurate stencil for the third derivative operator is correct and 
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
N = [100 200 400 800 1600 2000];
stencil = [-1, 3, -3 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute error for different dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err_d3f = zeros(size(N));
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
  d3f_exact = 27*sin(3*x).*exp(cos(3*x)) ...
            + 81*cos(3*x).*sin(3*x).*exp(cos(3*x)) ...
            - 27*(sin(3*x).^3).*exp(cos(3*x)) ...
            - 125*cos(5*x) + 336*x.^5;
  d3f_approx = zeros(size(x));
  for j = 4:length(x)-1
    d3f_approx(j) = sum(f(j-2:j+1).*stencil)/dx^3;
  end
  err_d3f(i) = norm(d3f_approx(4:end-3)-d3f_exact(4:end-3),'inf');

  % plot results
  figure(1); clf;
  plot(x,d3f_exact,'b');
  hold on;
  plot(x(4:end-3),d3f_approx(4:end-3),'rx');

end

% compute order of discretization
P = polyfit(log(N),log(err_d3f),1);
order = -P(1)

% plot error vs N
figure(2); clf;
loglog(N,err_d3f,'bo');
hold on;
N_plot = 100:1:10000;
plot(N_plot,exp(log(N_plot)*P(1)+P(2)),'r');
order_str = sprintf('Order = %1.1f', order);
text(1000,10,order_str);
