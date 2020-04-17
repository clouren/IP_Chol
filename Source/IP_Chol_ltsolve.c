//------------------------------------------------------------------------------
// IP_Chol/IP_Chol_ltsolve: Solve the system L' x = b for Cholesky
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE          \
return SLIP_OUT_OF_MEMORY;      \

#include "../Include/IP-Chol.h"

/* Purpose: This solves the system L'x = b for Cholesky factorization 
 * On input, L contains the lower triangular matrix. x has the solution
 * to the linear system from forward substitution
 */
SLIP_info IP_Chol_ltsolve 
(
    SLIP_matrix *L,     // The lower triangular matrix
    SLIP_matrix *x      // Solution vector
)
{
    SLIP_info ok;
    // Check input
    SLIP_REQUIRE(L, SLIP_CSC, SLIP_MPZ);
    SLIP_REQUIRE(x, SLIP_DENSE, SLIP_MPZ);
    
    int64_t p, j, n, k;
    // Set n
    n = L->n;
    // Iterate across the RHS vectors
    for (int64_t k = 0; k < x->n; k++)
    {
        // Iterate across the rows of x
        for (j = n-1; j >= 0; j--)
        {
            // if x(i,k) == 0 skip this operation
            if ( mpz_sgn ( SLIP_2D(x, j,k, mpz)) == 0) continue;
            
            for (p = L->p[j]+1; p < L->p[j+1]; p++)
            {
                OK ( SLIP_mpz_submul( SLIP_2D(x, j, k, mpz), L->x.mpz[p], 
                                      SLIP_2D( x, L->i[p], k, mpz)));
            }
            OK( SLIP_mpz_divexact( SLIP_2D(x, j, k, mpz), SLIP_2D(x, j, k, mpz), L->x.mpz[ L->p[j]]));
        }
    }
    return SLIP_OK;
}
