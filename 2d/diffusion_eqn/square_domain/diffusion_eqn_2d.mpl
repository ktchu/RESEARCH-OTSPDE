#
# This Maple script computes the source term and boundary conditions
# required in the diffusion equation 
#
#   u_t = kappa (u_xx + u_xx) + f
#
# on the domain [0,1]x[0,1] so that the exact solution is given by
# one of the following:
#
# (1)
# 
#   u(x,y,t) = 1 + y/3 + x/2 + xy/4 
#            - sin(pi*x)*sin(pi*y)*(1-exp(-2*kappa*pi^2*t)) 
#            + sin(2*pi*x)*sin(3*pi*y)*(1-exp(-13*kappa*pi^2*t)) 
#            + sin(5*pi*x)*sin(pi*y)*(2-exp(-26*kappa*pi^2*t)) 
#
# (2)
# 
#   u(x,y,t) = 1 + y/3 + x/2 + xy/4 
#            + sin(2*pi*x)*sin(3*pi*y)*exp(-5*kappa*pi^2*t) 
#            - sin(5*pi*x)*sin(pi*y)*exp(-2*kappa*pi^2*t) 
#
# Kevin T. Chu
# September 2007
#

#  Solution (1)
u(x,y,t) := 1 + y/3 + x/2 + x*y/4 
          - sin(pi*x)*sin(pi*y)*(1-exp(-2*kappa*pi^2*t))
          + sin(2*pi*x)*sin(3*pi*y)*(1-exp(-13*kappa*pi^2*t))
          + sin(5*pi*x)*sin(pi*y)*(2-exp(-26*kappa*pi^2*t));

u_t := diff(u(x,y,t),t);
u_xx := diff(diff(u(x,y,t),x),x);
u_yy := diff(diff(u(x,y,t),y),y);

f(x,y,t) := u_t - kappa*(u_xx+u_yy);
f(x,y,t) := simplify(f(x,y,t));

f_t(x,y,t) := diff(f(x,y,t),t);
f_t(x,y,t) := simplify(f_t(x,y,t));


#  Solution (2)
u(x,y,t) := 1 + y/3 + x/2 + x*y/4 
          + sin(2*pi*x)*sin(3*pi*y)*exp(-5*kappa*pi^2*t)
          - sin(5*pi*x)*sin(pi*y)*exp(-2*kappa*pi^2*t);

u_t := diff(u(x,y,t),t);
u_xx := diff(diff(u(x,y,t),x),x);
u_yy := diff(diff(u(x,y,t),y),y);

f(x,y,t) := u_t - kappa*(u_xx+u_yy);
f(x,y,t) := simplify(f(x,y,t));

f_t(x,y,t) := diff(f(x,y,t),t);
f_t(x,y,t) := simplify(f_t(x,y,t));

