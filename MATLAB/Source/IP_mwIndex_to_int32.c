//------------------------------------------------------------------------------
// IP_Chol/IP_mqIndex_to_int32: Convert an mwIndex to an int32_t
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"

/* Purpose: This function converts an mwIndex to an int32_t*/
void IP_mwIndex_to_int32
(
    int32_t* y,    // the int32_t vector
    mwIndex* x,    // the mwIndex vector
    mwSize n       // the size of x and y
)
{
    for (mwSize i = 0; i < n; i++)
    {  
        y[i] = (int32_t) x[i];
    }
}
