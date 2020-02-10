//------------------------------------------------------------------------------
// IP_Chol/IP_int32_to_mwIndex: Convert an int32_t to mwIndex
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"

/* Purpose: This function converts an int to an mwIndex */
void IP_int32_to_mwIndex
(
    mwIndex* y,    // the mwIndex vector
    int32_t* x,    // the int32_t vector
    int32_t n       // the size of x and y
)
{
    for (int32_t i = 0; i < n; i++)
    {        
        y[i] = (mwSize) x[i];
    }
}
