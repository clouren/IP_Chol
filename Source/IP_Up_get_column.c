//------------------------------------------------------------------------------
// IP_Chol/IP_Up_get_column: Obtain row k for up-looking factorization
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"


/* Purpose: This function obtains A(1:k-1,k) and stores it in a
 * dense vector x
 */
SLIP_info IP_Up_get_column  
(
    SLIP_sparse* A,    // input matrix
    int k,          // column to extract
    mpz_t* x        // A(1:k-1,k)
)
{
    SLIP_info ok;
    // Iterating accross the nonzeros in column k
    for (int i = A->p[k]; i < A->p[k+1]; i++)
    {
        if (A->i[i] <= k)
        {
            OK(SLIP_mpz_set(x[A->i[i]], A->x[i]));
        }
    }
    return ok;
}