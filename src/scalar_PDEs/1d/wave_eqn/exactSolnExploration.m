%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script was used to explore various possibilities for the choice of 
% exact solution for a 1D wave equation. 
%  
% Kevin T. Chu
% 2008 August
% Serendipity Research
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace
clear

% Physical parameters
c = 0.5;

% Cartesian grid
x_lo = -1.0;
x_hi =  1.0;
x = x_lo:(x_hi-x_lo)/100:x_hi;

% Generate movie of u as it varies over time
for t = 0:0.01:100

  % candidate #1
  % Initial conditions:
  %   u(x,0) = exp((-10*sin(0.5*pi*x))^2)
  %   u_t(x,0) = 0
  % Source term:
  %   f(x,t) = 0
%  u = 0.5*(exp(-10*(sin(0.5*pi*(x-c*t))).^2)+exp(-10*(sin(0.5*pi*(x+c*t))).^2));
%  plot_axes = [x_lo x_hi -2 2];

  % candidate #2
  % Initial conditions:
  %   u(x,0) = 2*exp((-10*sin(0.5*pi*x))^2)
  %   u_t(x,0) = -0.5*pi*sin(pi*x)
  % Source term:
  %   f(x,t) = 0
%  u = (exp(-10*(sin(0.5*pi*(x-c*t))).^2)+exp(-10*(sin(0.5*pi*(x+c*t))).^2)) ...
%    + 0.25*(cos(pi*(x+c*t)) - cos(pi*(x-c*t)));
%  plot_axes = [x_lo x_hi -3 3];

  % candidate #3
  % Initial conditions:
  %   u(x,0) = 2*exp((-10*sin(0.5*pi*x)^2)
  %   u_t(x,0) = -0.5*pi*sin(pi*x) + pi/c*cos(pi*x)
  % Source term:
  %   f(x,t) = 0.125/c*pi^2*sin(pi*(x-t/4))
  u = (exp(-10*(sin(0.5*pi*(x-c*t))).^2)+exp(-10*(sin(0.5*pi*(x+c*t))).^2)) ...
    + 0.25/c*(cos(pi*(x+c*t)) - cos(pi*(x-c*t))) ...
    + -0.25/c/(c+0.25)*( sin(pi*(x-0.25*t)) - sin(pi*(x+c*t)) ) ...
    +  0.25/c/(c-0.25)*( sin(pi*(x-0.25*t)) - sin(pi*(x-c*t)) );
  plot_axes = [x_lo x_hi -5 5];

  plot(x,u);
  axis(plot_axes);
  drawnow;
end
