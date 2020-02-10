//------------------------------------------------------------------------------
// IP_Chol/IP_mex_check_for_inf: Check A and b for infinity
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"


/* Purpose: This function checks if the input matrix or RHS has numbers too
 * large for double*/
void IP_mex_check_for_inf
(
    double* x, // The array of numeric values 
    mwSize n   // size of array
)
{
    for (mwSize k = 0; k < n; k++)
    {
        if (mxIsInf(x[k]))
        {
            mexErrMsgTxt("Numbers are too large for double. Please try the C "
                "code with mpfr, mpq, or mpz");
        }
    }
}
