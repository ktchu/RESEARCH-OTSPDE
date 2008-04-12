# 
# This Maple script computes the stencil for the second-order and
# fourth-order accurate central discretizations of the third-derivative
# operator in one space dimension.
#
# NOTE:
# - The linear algebra package MUST be loaded before this  
#   computation.  This can be done using the 'with(linalg)' 
#   command.
# 
# Kevin T. Chu 
# September 2007
# 

# construct coefficient matrix for 2nd-order accurate stencil
# (from Taylor series expansion)
C_2 := matrix([ [4   , 2  ], 
                [8/3 , 1/3] ]);

# construct RHS that picks out the third-derivative term
RHS_2 := vector([0, 1]);

# compute stencil
stencil := linsolve(C_2,RHS_2);

# compute coefficient on O(h^5) term
error_coef_5 := 2/5!*(32*stencil[1]+stencil[2]);


# construct coefficient matrix for 4th-order accurate stencil
# (from Taylor series expansion)
# NOT DONE YET
