#
# This Maple script computes the derivative of the cubic interpolant
# used to locate the interface defined by the zero level set of phi.
# The values of phi are given at adjacent grid points: x1 < x2 < x3 < x4.
# phi changes signs between x2 and x3.  The spacing between grid points 
# dx.
#
# Kevin T. Chu
# September 2007
#

phi(x) := -1/6*phi_1*(x-x2)*(x-x3)*(x-x4)
        + 1/2*phi_2*(x-x1)*(x-x3)*(x-x4)
        - 1/2*phi_3*(x-x1)*(x-x2)*(x-x4)
        + 1/6*phi_4*(x-x1)*(x-x2)*(x-x3);

dphidx := diff(phi(x),x);
