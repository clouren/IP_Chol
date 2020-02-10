//------------------------------------------------------------------------------
// IP_Chol/IP_gmp_mex_realloc: A gmp realloc function for matlab
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"

/* Purpose: A GMP reallocation function 
 * This allows GMP to use MATLAB's default realloc function 
 */

// A GMP realloc function
void* IP_gmp_mex_realloc 
(
    void* x,    // void* to be reallocated 
    size_t a,   // Previous size
    size_t b    // New size
)
{
    return mxRealloc(x,b);
}
