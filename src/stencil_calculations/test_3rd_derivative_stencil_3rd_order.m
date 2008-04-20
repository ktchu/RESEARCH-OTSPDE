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
N = [25 50 100 200 400 800];
stencil = [-7/4, 41/4, -49/2, 59/2, -71/4, 17/4];  % fully upwind (-5:0)
stencil = [1/4, -7/4, 7/2, -5/2, 1/4, 1/4];  % partial upwind (-3:2)
stencil = [-1/4, -1/4, 5/2, -7/2, 7/4, -1/4];  % partial upwind (-2:3)
stencil = [-7/4, 25/4, -17/2, 11/2, -7/4, 1/4]; % partial downwind (-1:4)


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

%  for partial upwind (-3:2)
%  for j = 4:length(x)-3
%    d3f_approx(j) = sum(f(j-3:j+2).*stencil)/dx^3;
%  end
%  err_d3f(i) = norm(d3f_approx(4:end-3)-d3f_exact(4:end-3),'inf');

%  for partial upwind (-2:3)
%  for j = 3:length(x)-3
%    d3f_approx(j) = sum(f(j-2:j+3).*stencil)/dx^3;
%  end
%  err_d3f(i) = norm(d3f_approx(3:end-3)-d3f_exact(3:end-3),'inf');

%  for partial downwind(-1:4)
  for j = 2:length(x)-4
    d3f_approx(j) = sum(f(j-1:j+4).*stencil)/dx^3;
  end
  err_d3f(i) = norm(d3f_approx(2:end-4)-d3f_exact(2:end-4),'inf');

%  for fully upwind (-5:0)
%  for j = 6:length(x)
%    d3f_approx(j) = sum(f(j-5:j).*stencil)/dx^3;
%  end
%  err_d3f(i) = norm(d3f_approx(6:end)-d3f_exact(6:end),'inf');

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
N_plot = 10:1:1000;
plot(N_plot,exp(log(N_plot)*P(1)+P(2)),'r');
order_str = sprintf('Order = %1.1f', order);
text(100,10,order_str);
