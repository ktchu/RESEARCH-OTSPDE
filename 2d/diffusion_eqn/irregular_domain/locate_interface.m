%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB function locates the position of the zero level set
% of a cubic interpolant of a function phi along one direction given
% the value of phi at four positions and an initial guess for the 
% position of the zero level set.
%
% A simple Newton iteration is used to refine the initial guess.
%
% Kevin T. Chu
% 2007 September
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x0 = locate_interface(phi,X,x0)

dx = abs(X(2)-X(1));
max_iters = 20;
tol = 1e-10;

x_start = x0;

count = 1;
diffs = x0 - X;
res = 1/dx^3*(-1/6*phi(1)*diffs(2)*diffs(3)*diffs(4) ...
             + 1/2*phi(2)*diffs(1)*diffs(3)*diffs(4) ...
             - 1/2*phi(3)*diffs(1)*diffs(2)*diffs(4) ...
             + 1/6*phi(4)*diffs(1)*diffs(2)*diffs(3));
while (abs(res) > tol & count < max_iters)

  diffs = x0 - X;
  res = 1/dx^3*(-1/6*phi(1)*diffs(2)*diffs(3)*diffs(4) ...
               + 1/2*phi(2)*diffs(1)*diffs(3)*diffs(4) ...
               - 1/2*phi(3)*diffs(1)*diffs(2)*diffs(4) ...
               + 1/6*phi(4)*diffs(1)*diffs(2)*diffs(3));

  J = 1/dx^3 ...
     *(-1/6*phi(1)*(diffs(2)*diffs(3)+diffs(2)*diffs(4)+diffs(3)*diffs(4)) ...
      + 1/2*phi(2)*(diffs(1)*diffs(3)+diffs(1)*diffs(4)+diffs(3)*diffs(3)) ...
      - 1/2*phi(3)*(diffs(1)*diffs(2)+diffs(1)*diffs(4)+diffs(2)*diffs(4)) ...
      + 1/6*phi(4)*(diffs(1)*diffs(2)+diffs(1)*diffs(3)+diffs(2)*diffs(3)));

  x0 = x0 - res/J;

  count = count + 1;

end
