//------------------------------------------------------------------------------
// IP_Chol/IP_array_div: divide a vector by a scalar
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function takes as input a mpz_t** array and divides it by a
 * mpz_t constant storing the solution in a mpq_t** array. This is used 
 * internally to divide the solution vector by the determinant of the matrix.
 * 
 * On output, the contents of the array x2 are modified
 * 
 */

#define FREE_WORKSPACE       \
    SLIP_MPQ_CLEAR(det2);    \

# include "../Include/IP-Chol.h"

SLIP_info IP_array_div // divides the x vector by the determinant
(
    mpq_t** x2,     // solution of x/det
    mpz_t** x,      // input vector
    mpz_t det,      // given determinant of matrix
    int32_t n,      // size of x and x2 
    int32_t numRHS  // number of rhs vectors
)
{
    SLIP_info ok;
    if (!x2 || !x) {return SLIP_INCORRECT_INPUT;}
    // Set det2 = det
    mpq_t det2; SLIP_MPQ_SET_NULL(det2);
    OK(SLIP_mpq_init(det2)); 
    OK(SLIP_mpq_set_num(det2, det));

    //--------------------------------------------------------------------------
    // iterate each entry of x, copy to x2 and divide it by det
    //--------------------------------------------------------------------------
    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t k = 0; k < numRHS; k++)
        {
            // Set x2[i] = x[i]
            OK(SLIP_mpq_set_num(x2[i][k], x[i][k]));
            // x2[i] = x2[i] / det2
            OK(SLIP_mpq_div(x2[i][k], x2[i][k], det2));
        }
    }

    //--------------------------------------------------------------------------
    // Free memory associated with det2
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
    return SLIP_OK;
}

