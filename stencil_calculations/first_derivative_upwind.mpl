# 
# This Maple script computes the stencil for the first-, third-, and
# sixth-order accurate upwind discretizations of the first-derivative
# operator in one space dimension when the upwind direction is in the
# negative x-direction.
#
# NOTE:
# - The linear algebra package MUST be loaded before this  
#   computation.  This can be done using the 'with(linalg)' 
#   command.
# 
# Kevin T. Chu 
# September 2007
# 

with(linalg);

# construct coefficient matrix for 1st-order accurate stencil
# (from Taylor series expansion)
C_1 := matrix([ [1 , 1],
                [0 , -1] ]);

# construct RHS that picks out the first-derivative term
RHS_1 := vector([0, 1]);

# compute stencil
stencil_1 := linsolve(C_1,RHS_1);

# compute coefficient on O(h^2) term
error_coef_2 := 1/2!*stencil_1[2];

# construct coefficient matrix for 3rd-order accurate stencil
# (from Taylor series expansion)
C_3 := matrix([ [1     , 1 , 1      , 1    ],
                [1     , 0 , -1     , -2   ],
                [1/2   , 0 , 1/2    , 2    ],
                [1/6   , 0 , -1/6   , -4/3 ] ]);

# construct RHS that picks out the first-derivative term
RHS_3 := vector([0, 1, 0, 0]);

# compute stencil
stencil_3 := linsolve(C_3,RHS_3);

# compute coefficient on O(h^4), O(h^5), O(h^6) terms
error_coef_4 := 1/4!*(stencil_3[1]+stencil_3[3]+16*stencil_3[4]);
error_coef_5 := 1/5!*(stencil_3[1]-stencil_3[3]-32*stencil_3[4]);
error_coef_6 := 1/6!*(stencil_3[1]+stencil_3[3]+64*stencil_3[4]);

# construct coefficient matrix for 5th-order accurate stencil
# (from Taylor series expansion) 
C_5 := matrix([ [1    , 1     , 1 , 1      , 1     , 1     ],
                [2    , 1     , 0 , -1     , -2    , -3    ],
                [2    , 1/2   , 0 , 1/2    , 2     , 9/2   ],
                [4/3  , 1/6   , 0 , -1/6   , -4/3  , -9/2  ],
                [2/3  , 1/24  , 0 , 1/24   , 2/3   , 27/8  ],
                [4/15 , 1/120 , 0 , -1/120 , -4/15 , -81/40] ]);

# construct RHS that picks out the first-derivative term
RHS_5 := vector([0, 1, 0, 0, 0, 0]);

# compute stencil
stencil_5 := linsolve(C_5,RHS_5);

# compute coefficient on O(h^6) term
error_coef_6 := 1/6!*(3^6*stencil_5[1]+2^6*stencil_5[2]+stencil_5[3]
                     +stencil_5[5]+2^6*stencil_5[6]);


# construct coefficient matrix for 6th-order accurate stencil
# (from Taylor series expansion ... only odd terms are kept to reduce
#  size of matrix)
C_6 := matrix([ [6     , 4    , 2    ],
                [9     , 8/3  , 1/3  ],
                [81/20 , 8/15 , 1/60 ] ]);

# construct RHS that picks out the first-derivative term
RHS_6 := vector([1, 0, 0]);

# compute stencil
stencil_6 := linsolve(C_6,RHS_6);

# compute coefficient on O(h^7) term
error_coef_7 := 2/7!*(3^7*stencil_6[1]+2^7*stencil_6[2]+stencil_6[3]);

