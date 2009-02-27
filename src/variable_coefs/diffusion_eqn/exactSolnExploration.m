%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script was used to explore various possibilities for the choice of 
% exact solution for a 1D diffusion equation with a spatially varying
% diffusivity.
%
%   u_t = d/dx ( D(x) * u_x ) + f(x,t)
%  
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
  %   u(x,0) = ??
  % Source term:
  %   f(x,t) = ??
  u = exp(sin(2*pi*(x-t))) - exp(4*cos(pi*(x-t)));
  plot_axes = [x_lo x_hi -20 5];

  % candidate #2
  % Initial conditions:
  %   u(x,0) = ??
  % Source term:
  %   f(x,t) = ??
  u = exp(sin(2*pi*(x-t)) + cos(3*pi*(x-t)));
  plot_axes = [x_lo x_hi 0 10];

  % candidate #3
  % Diffusivity:
  %   D(x) = 1 + 0.5*cos(3*pi*x)
  % Initial conditions:
  %   u(x,0) = exp(sin(2*pi*x) + cos(3*pi*x))
  % Source term:
  %   f(x,t) = 0.5*pi*exp(sin(2*pi*(x-t))+cos(3*pi*x)) ...
  %          .*( (-27*pi + 9*pi*cos(3*pi*x) + 36*pi*cos(3*pi*x).^2 ...
  %              + 9*pi*cos(3*pi*x).^3) ...
  %            + cos(2*pi*(x-t)).*(-4 + 30*pi*sin(3*pi*x) ...
  %                               + 12*pi*cos(3*pi*x).*sin(3*pi*x)) ...
  %            + (sin(2*pi*(x-t)) - cos(2*pi*(x-t)).^2) ...
  %            .*(8*pi + 4*pi*cos(3*pi*x)) );
  u = exp(sin(2*pi*(x-t)) + cos(3*pi*x));
  plot_axes = [x_lo x_hi 0 10];

  % candidate #4
  % Diffusivity:
  %   D(x) = 2-x.^2
  % Initial conditions:
  %   u(x,0) = (1-x.^2).*(1+0.5*sin(3*pi*x))
  % Source term:
  %   f(x,t) = (4-6*x.^2) ...
  %          + cos(3*pi*(x-t)).*( 1.5*pi*(x.^2-1) + 3*pi*(5*x-3*x.^3) ) ...
  %          + sin(3*pi*(x-t)).*( 2-3*x.^2+9*pi^2*(1-1.5*x.^2+0.5*x.^4) )
  u = (1-x.^2).*(1+0.5*sin(3*pi*(x-t)));
  plot_axes = [x_lo x_hi 0 2];

  plot(x,u);
  axis(plot_axes);
  drawnow;
end
