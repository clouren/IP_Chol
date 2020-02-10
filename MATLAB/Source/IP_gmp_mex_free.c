//------------------------------------------------------------------------------
// IP_Chol/IP_gmp_mex_free: A gmp free function for matlab
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"

/* Purpose: A GMP free function. This allows GMP to use
 * MATLAB's mxFree instead of free */

// A GMP realloc function
void IP_gmp_mex_free
(
    void* x,    // void* to be freed
    size_t a    // Size
)
{
    if (x) 
        mxFree(x);
}
