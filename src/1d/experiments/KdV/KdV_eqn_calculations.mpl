# 
# This Maple script computes u_tt for the KdV equations.
#
# Kevin T. Chu 
# September 2007
# 

# u_x
u_x(x) := diff(u(x),x);

# u_t
u_t(x) := -diff(diff(diff(u(x),x),x),x) - 6*u(x)*u_x(x);

# u_tx
u_tx(x) := diff(u_t(x),x);

# u_txxx
u_txxx(x) := diff(diff(diff(u_t(x),x),x),x);

# u_tt
u_tt(x) := -u_txxx(x) - 6*u_t(x)*u_x(x) - 6*u(x)*u_tx(x);
u_tt(x) := simplify(u_tt(x));

# O(dt^2) term for semi-implicit scheme
correction_semi_implicit := u_txxx(x) - 6*u_t(x)*u_x(x) - 6*u(x)*u_tx(x);
correction_semi_implicit := simplify(correction_semi_implicit);
