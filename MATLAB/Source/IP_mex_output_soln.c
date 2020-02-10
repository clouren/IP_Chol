//------------------------------------------------------------------------------
// IP_Chol/IP_mex_output_soln: Output x to matlab
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"


/* Purpose: This function outputs the solution as a mxArray. */
mxArray* IP_mex_output_soln
(
    double** x,        // The matrix to be output
    int32_t m,         // size of x
    int32_t n          // the size of x 
)
{
    // Create a m * n array
    mxArray* Xmatlab = mxCreateDoubleMatrix ((mwSize) m, (mwSize) n, mxREAL);

    // Get the numeric values
    double* x2 = mxGetPr(Xmatlab);
    int32_t count = 0;

    // Populate the nonzeros in output matrix
    for (int32_t j = 0; j < n; j++)
    {
        for (int32_t i = 0; i < m; i++)
        {
            x2[count] = x[i][j];
            count+=1;
        }
    }
    return Xmatlab;
}
