%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script was used to explore various possibilities for the choice of 
% exact solution for a 1D wave equation with a variable wave speed.
%
%   u_tt - (c(x))^2 * u_xx = f(x,t)
%  
% with c(x) = 1 + 0.5*sin(pi*x/2).
% 
% Kevin T. Chu
% 2008 August
% Serendipity Research
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear workspace
clear

% Physical parameters

% Cartesian grid
x_lo = -1.0;
x_hi =  1.0;
x = x_lo:(x_hi-x_lo)/100:x_hi;

% Generate movie of u as it varies over time
for t = 0:0.01:100

  % candidate #1
  % Initial conditions:
  %   u(x,0) = 2*exp((-10*sin(0.5*pi*x))^2) + 0.5*sin(pi*x) - 0.5*cos(pi*x)
  %   u_t(x,0) = -0.5*pi*(sin(pi*x) - cos(pi*x))
  % Source term:
  %   f(x,t) = pi^2/2*cos(3*pi*x) ...
  %    * (-2*cos(pi*(x-t) + 2*sin(pi*(x+t)) ...
  %      - 180*cos(3*pi*x)*cos(0.5*pi*(x-t))^2*exp(-10*sin(0.5*pi*(x-t))^2) ...
  %      - 180*cos(3*pi*x)*cos(0.5*pi*(x+t))^2*exp(-10*sin(0.5*pi*(x+t))^2) ...
  %      + 200*cos(3*pi*x)*cos(0.5*pi*(x-t))^4*exp(-10*sin(0.5*pi*(x-t))^2) ...
  %      + 200*cos(3*pi*x)*cos(0.5*pi*(x+t))^4*exp(-10*sin(0.5*pi*(x+t))^2) ...
  %      - 360*cos(0.5*pi*(x-t))^2*exp(-10*sin(0.5*pi*(x-t))^2) ...
  %      - 360*cos(0.5*pi*(x+t))^2*exp(-10*sin(0.5*pi*(x+t))^2) ...
  %      + 400*cos(0.5*pi*(x-t))^4*exp(-10*sin(0.5*pi*(x-t))^2) ...
  %      + 400*cos(0.5*pi*(x+t))^4*exp(-10*sin(0.5*pi*(x+t))^2) ...
  %      + cos(3*pi*)*sin(pi*(x+t)) - cos(3*pi*)*cos(pi*(x-t)) ...
  %      - 20*exp(-10*sin(0.5*pi*(x-t))^2) ...
  %      - 20*exp(-10*sin(0.5*pi*(x+t))^2) ...
  %      - 10*cos(3*pi*x)*exp(-10*sin(0.5*pi*(x-t))^2) ...
  %      - 10*cos(3*pi*x)*exp(-10*sin(0.5*pi*(x+t))^2); 
  %
  u = (exp(-10*(sin(0.5*pi*(x-t))).^2) + exp(-10*(sin(0.5*pi*(x+t))).^2)) ...
    + 0.5*(sin(pi*(x+t)) - cos(pi*(x-t)));
  plot_axes = [x_lo x_hi -3 3];

  % candidate #2
  % Initial conditions:
  %   u(x,0) = 3*exp((-10*sin(0.5*pi*x))^2) 
  %   u_t(x,0) = -10*pi*sin(0.5*pi*x)*cos(0.5*pi*x) ...
  %            * exp(-10*sin(0.5*pi*x)^2)
  % Source term:
  %   f(x,t) = ??
  u = (exp(-10*(sin(0.5*pi*(x-t))).^2) + 2*exp(-10*(sin(0.5*pi*(x+t))).^2));
  plot_axes = [x_lo x_hi -3 3];

  % candidate #3
  % Initial conditions:
  %   u(x,0) = sin(2*pi*x) + cos(3*pi*x)
  %   u_t(x,0) = -pi*(2*cos(2*pi*x) + 3*sin(3*pi*x))
  % Source term:
  %   f(x,t) = 0.25*pi^2*sin(0.5*pi*x) ...
  %          * ( 16*sin(2*pi*(x-t)) + 36*cos(3*pi*(x+t)) ...
  %            + 4*sin(0.5*pi*x)*sin(2*pi*(x-t)) ...
  %            + 9*sin(0.5*pi*x)*cos(3*pi*(x+t)) )
  u = sin(2*pi*(x-t)) + cos(3*pi*(x+t));
  plot_axes = [x_lo x_hi -3 3];

  plot(x,u);
  axis(plot_axes);
  drawnow;
end
