#
# This Maple script computes the source function required in order to 
# have to solution to the fourth-order parabolic equation:
#
#   u_t = -(u_xxxx + 2 u_xxyy + u_yyyy) + f(x,y,t)
#
# be given by
#
#   u(x,y,t) = sin(2*pi*x)*sin(5*pi*y)*cos(3*t)
#            + sin(pi*x)*sin(3*pi*y)*(1-exp(-10*pi^2*t))
#
# Kevin T. Chu
# September 2007
#

u(x,y,t) := sin(2*pi*x)*sin(5*pi*y)*cos(3*t)
          + sin(pi*x)*sin(3*pi*y)*(1-exp(-10*pi^2*t));

u_t(x,y,t) := diff(u(x,y,t),t);
u_xxxx(x,y,t) := diff(diff(diff(diff(u(x,y,t),x),x),x),x);
u_xxyy(x,y,t) := diff(diff(diff(diff(u(x,y,t),x),x),y),y);
u_yyyy(x,y,t) := diff(diff(diff(diff(u(x,y,t),y),y),y),y);

f(x,y,t) := u_t(x,y,t) + (u_xxxx(x,y,t) + 2*u_xxyy(x,y,t) + u_yyyy(x,y,t));
