//------------------------------------------------------------------------------
// IP_Chol/IP_array_mul: multiply a vector by a scalar
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* 
 * Purpose: This function multiplies vector x by the determinant of matrix. 
 * 
 * On output the contents of the x vector is modified 
 */
SLIP_info IP_array_mul // multiplies vector x by the determinant of matrix
(
    mpz_t** x,      // matrix to be multiplied
    mpz_t det,      // given determinant of matrix
    int32_t n,      // size of x
    int32_t numRHS  // number of RHS vectors
)
{
    SLIP_info ok;
    if (!x) {return SLIP_INCORRECT_INPUT;}
    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t k = 0; k < numRHS; k++)
        {
            // x[i][k] = x[i][k] * det
            OK(SLIP_mpz_mul(x[i][k], x[i][k], det));
        }
    }
    return SLIP_OK;
}

