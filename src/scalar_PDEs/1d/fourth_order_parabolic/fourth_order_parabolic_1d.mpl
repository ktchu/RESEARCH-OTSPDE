#
# This Maple script computes the source function required in order to 
# have to solution to the fourth-order parabolic equation:
#
#   u_t = -u_xxxx + f(x,y,t)
#
# be given by
#
#   u(x,t) = 2*sin(2*Pi*x)*cos(10*Pi*t)
#          + sin(5*Pi*x)*(2-exp(-10*Pi^2*t))
#
# Kevin T. Chu
# October 2007
#

u(x,t) := 2*sin(2*Pi*x)*cos(10*Pi*t)
        + sin(5*Pi*x)*(2-exp(-10*Pi^2*t));

u_t(x,t) := diff(u(x,t),t);
u_xxxx(x,t) := diff(diff(diff(diff(u(x,t),x),x),x),x);

f(x,t) := u_t(x,t) + u_xxxx(x,t);
f_t(x,t) := simplify(diff(f(x,t),t));
