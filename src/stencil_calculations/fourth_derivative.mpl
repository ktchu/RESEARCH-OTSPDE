# 
# This Maple script computes the stencil for the second- and fourth-order 
# accurate discretizations of the fourth-derivative operator in one space 
# dimension.
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
C_2 := matrix([ [2    , 2     , 1],  
                [4    , 1     , 0],
                [4/3  , 1/12  , 0] ]);

# construct RHS that picks out the fourth-derivative term
RHS_2 := vector([0, 0, 1]);

# compute stencil
stencil := linsolve(C_2,RHS_2);

# compute coefficient on O(h^6) term
error_coef_6 := 2/6!*(64*stencil[1]+stencil[2]);

# construct coefficient matrix (from Taylor series expansion)
C_4 := matrix([ [2    , 2    , 2     , 1],  
                [9    , 4    , 1     , 0],
                [27/4 , 4/3  , 1/12  , 0],
                [81/40  , 8/45 , 1/360 , 0] ]);

# construct RHS that picks out the fourth-derivative term
RHS_4 := vector([0, 0, 1, 0]);

# compute stencil
stencil := linsolve(C_4,RHS_4);

# compute coefficient on leading order error
error_coef := 2/8!*(3^8*stencil[1]+256*stencil[2]+stencil[3]);

