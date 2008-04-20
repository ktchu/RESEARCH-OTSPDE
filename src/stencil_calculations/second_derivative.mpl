# 
# This Maple script computes the stencil for the second- and fourth-order 
# accurate accurate upwind discretizations of the second-derivative operator 
# in one space dimension. 
#
# NOTE:
# - The linear algebra package MUST be loaded before this  
#   computation.  This can be done using the 'with(linalg)' 
#   command.
# 
# Kevin T. Chu 
# September 2007
# 

# construct coefficient matrix for second-order accurate stencil
# (from Taylor series expansion)
C_2 := matrix([ [1     , 1 , 1      ],
                [1     , 0 , -1     ],
                [1/2   , 0 , 1/2    ] ] );

# construct RHS that picks out the second-derivative term
RHS_2 := vector([0, 0, 1]);

# compute stencil
stencil_2 := linsolve(C_2,RHS_2);

# compute coefficient on O(h^4) term
error_coef_4 := 1/4!*(stencil_2[1]+stencil_2[3]);

# construct coefficient matrix for fourth-order accurate stencil
# (from Taylor series expansion)
C_4 := matrix([ [1    , 1     , 1 , 1      , 1     ],
                [2    , 1     , 0 , -1     , -2    ],
                [2    , 1/2   , 0 , 1/2    , 2     ],
                [4/3  , 1/6   , 0 , -1/6   , -4/3  ],
                [2/3  , 1/24  , 0 , 1/24   , 2/3   ],
                [4/15 , 1/120 , 0 , -1/120 , -4/15 ] ]);

# construct RHS that picks out the second-derivative term
RHS_4 := vector([0, 0, 1, 0, 0, 0]);

# compute stencil
stencil_4 := linsolve(C_4,RHS_4);

# compute coefficient on O(h^6) term
error_coef_6 := 1/6!*(64*stencil_4[1]+stencil_4[2]+stencil_4[4]
                     +64*stencil_4[5]);

