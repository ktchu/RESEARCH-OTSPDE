# 
# This Maple script computes the stencil for the first-order and
# third-order accurate upwind discretizations of the third-derivative
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

#################################################################
#
# This section computes partial upwind stencils (-3:2)
#
#################################################################

with(linalg);

# construct coefficient matrix for 1st-order accurate stencil
# (from Taylor series expansion)
C_1 := matrix([ [1   , 1 , 1    , 1],
                [1   , 0 , -1   , -2],
                [1/2 , 0 , 1/2  , 2],
                [1/6 , 0 , -1/6 , -4/3] ]);

# construct RHS that picks out the third-derivative term
RHS_1 := vector([0, 0, 0, 1]);

# compute stencil
stencil_1 := linsolve(C_1,RHS_1);

# compute coefficient on O(h^4) term
error_coef_4 := 1/4!*(stencil_1[1]+stencil_1[3]+2/3*stencil_1[4]);

# construct coefficient matrix for 3rd-order accurate stencil
# (from Taylor series expansion)
C_3 := matrix([ [1    , 1     , 1 , 1      , 1     , 1     ],
                [2    , 1     , 0 , -1     , -2    , -3    ],
                [2    , 1/2   , 0 , 1/2    , 2     , 9/2   ],
                [4/3  , 1/6   , 0 , -1/6   , -4/3  , -9/2  ],
                [2/3  , 1/24  , 0 , 1/24   , 2/3   , 27/8  ],
                [4/15 , 1/120 , 0 , -1/120 , -4/15 , -81/40] ]);

# construct RHS that picks out the third-derivative term
RHS_3 := vector([0, 0, 0, 1, 0, 0]);

# compute stencil
stencil_3 := linsolve(C_3,RHS_3);

# compute coefficient on O(h^6) term
error_coef_6 := 1/6!*(64*stencil_3[1]+stencil_3[2]+stencil_3[4]
                     +64*stencil_3[5]+729*stencil_3[6]);


#################################################################
#
# This section computes partial upwind stencils (-2:3)
#
#################################################################

# construct coefficient matrix for 3rd-order accurate stencil
# (from Taylor series expansion)
C_3 := matrix([ [1     , 1    , 1     , 1 , 1      , 1     ],
                [3     , 2    , 1     , 0 , -1     , -2    ],
                [9/2   , 2    , 1/2   , 0 , 1/2    , 2     ],
                [9/2   , 4/3  , 1/6   , 0 , -1/6   , -4/3  ],
                [27/8  , 2/3  , 1/24  , 0 , 1/24   , 2/3   ],
                [81/40 , 4/15 , 1/120 , 0 , -1/120 , -4/15 ] ]);

# construct RHS that picks out the third-derivative term
RHS_3 := vector([0, 0, 0, 1, 0, 0]);

# compute stencil
stencil_3 := linsolve(C_3,RHS_3);

# compute coefficient on O(h^6) term
error_coef_6 := 1/6!*(729*stencil_3[1]+64*stencil_3[2]+stencil_3[3]
                     +stencil_3[5]+64*stencil_3[6]);


#################################################################
#
# This section computes partial upwind stencils (-1:4)
#
#################################################################

# construct coefficient matrix for 3rd-order accurate stencil
# (from Taylor series expansion)
C_3 := matrix([ [1      , 1    , 1     , 1     , 1 , 1     ],
                [4      , 3     , 2    , 1     , 0 , -1    ],
                [8      , 9/2   , 2    , 1/2   , 0 , 1/2   ],
                [32/3   , 9/2   , 4/3  , 1/6   , 0 , -1/6  ],
                [32/3   , 27/8  , 2/3  , 1/24  , 0 , 1/24  ],
                [128/15 , 81/40 , 4/15 , 1/120 , 0 , -1/120] ]);

# construct RHS that picks out the third-derivative term
RHS_3 := vector([0, 0, 0, 1, 0, 0]);

# compute stencil
stencil_3 := linsolve(C_3,RHS_3);

# compute coefficient on O(h^6) term
error_coef_6 := 1/6!*(4096*stencil_3[1]+729*stencil_3[2]+64*stencil_3[3]
                     +stencil_3[4]+stencil_3[6]);


#################################################################
#
# This section computes fully upwind stencils (-5:0)
#
#################################################################

# construct coefficient matrix for 1st-order accurate stencil
# (from Taylor series expansion)
C_1 := matrix([ [1 , 1    , 1    , 1   ],
                [0 , -1   , -2   , -3  ],
                [0 , 1/2  , 2    , 9/2 ],
                [0 , -1/6 , -4/3 , -9/2] ]);

# construct RHS that picks out the third-derivative term
RHS_1 := vector([0, 0, 0, 1]);

# compute stencil
stencil_1 := linsolve(C_1,RHS_1);

# compute coefficient on O(h^4) term
error_coef_4 := 1/4!*(stencil_1[2]+16*stencil_1[3]+81*stencil_1[4]);

# construct coefficient matrix for 3rd-order accurate stencil
# (from Taylor series expansion)
C_3 := matrix([ [1 , 1      , 1     , 1      , 1       , 1       ],
                [0 , -1     , -2    , -3     , -4      , -5      ],
                [0 , 1/2    , 2     , 9/2    , 8       , 25/2    ],
                [0 , -1/6   , -4/3  , -9/2   , -32/3   , -125/6  ],
                [0 , 1/24   , 2/3   , 27/8   , 32/3    , 625/24  ],
                [0 , -1/120 , -4/15 , -81/40 , -128/15 , -625/24 ] ]);

# construct RHS that picks out the third-derivative term
RHS_3 := vector([0, 0, 0, 1, 0, 0]);

# compute stencil
stencil_3 := linsolve(C_3,RHS_3);

# compute coefficient on O(h^6) term
error_coef_6 := 1/6!*(stencil_3[2]+2^6*stencil_3[3]+3^6*stencil_3[4]
                     +4^6*stencil_3[5]+5^6*stencil_3[6]);

