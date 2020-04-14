//------------------------------------------------------------------------------
// IP_Chol/IP_Chol_ltsolve: Solve the system L' x = b for Cholesky
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE          \
return SLIP_OUT_OF_MEMORY;      \

#include "../Include/IP-Chol.h"

/* Purpose: This solves the system L'x = b for Cholesky factorization */
int IP_Chol_ltsolve 
(
    SLIP_matrix *L,    // The lower triangular matrix
    SLIP_matrix *x      // Solution vector
)
{
    SLIP_info ok;
    int p, j, n, k;
    n = L->n;
    for (int k = 0; k < numRHS; k++)
    {
        for (j = n-1; j >= 0; j--)
        {
            if ( mpz_sgn ( SLIP_2D(x, i,j, mpz)) == 0) continue;
            
            for (p = L->p[j]+1; p < L->p[j+1]; p++)
            {
                OK ( SLIP_mpz_submul( SLIP_2D(x, j, k, mpz), L->x[p], 
                                      SLIP_2D( x, L->i[p], k, mpz)));
            }
            OK( SLIP_mpz_divexact( SLIP_2D(x, j, k, mpz), SLIP_2D(x, j, k, mpz), L->x[ L->p[j]]));
        }
    }
    return SLIP_OK;
}
