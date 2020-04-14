///------------------------------------------------------------------------------
// IP_Chol/IP_forward_sub: Solve the system LDx = b
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function performs sparse REF forward substitution This is
 * essentially the same as the sparse REF triangular solve applied to each
 * column of the right hand side vectors. Like the normal one, this
 * function expects that the vector x is dense. As a result,the nonzero
 * pattern is not computed and each nonzero in x is iterated across.
 * The system to solve is LDx = x
 *
 * On output, the mpz_t** x structure is modified
 *
 */

#define SLIP_FREE_WORKSPACE            \
{                                      \
    for (i = 0; i < n; i++)            \
    {                                  \
        SLIP_FREE(h[i]);               \
    }                                  \
    SLIP_FREE(h);                      \
}

#include "../Include/IP-Chol.h"

SLIP_info IP_forward_sub
(
    SLIP_matrix *L,   // lower triangular matrix
    SLIP_matrix *x,        // right hand side matrix of size n*numRHS
    SLIP_matrix *rhos,      // sequence of pivots used in factorization
    int32_t numRHS    // number of columns in x
)
{
    SLIP_info ok;
    int32_t i, j, p, k, n, m, mnew, sgn, **h;
    // Size of x vector
    n = L->n;

    // calloc is used, so that h is initialized for SLIP_FREE_WORKSPACE
    h = (int32_t**) SLIP_calloc(n, sizeof(int32_t*));
    if (!h)
    {
        return SLIP_OUT_OF_MEMORY;
    }
    for (i = 0; i < n; i++)
    {
        h[i] = (int32_t*) SLIP_malloc(numRHS* sizeof(int32_t));
        if (!h[i])
        {
            SLIP_FREE_WORKSPACE;
            return SLIP_OUT_OF_MEMORY;
        }
        for (j = 0; j < numRHS; j++)
        {
            h[i][j] = -1;
        }
    }

    //--------------------------------------------------------------------------
    // Iterate across each RHS vector
    //--------------------------------------------------------------------------

    for (k = 0; k < numRHS; k++)
    {
        //----------------------------------------------------------------------
        // Iterate accross all nonzeros in x. Assume x is dense
        //----------------------------------------------------------------------
        for (i = 0; i < n; i++)
        {
            p = h[i][k];
            // If x[i][k] = 0, can skip operations and continue to next i
            OK(SLIP_mpz_sgn(&sgn, SLIP_2D(x, i, k, mpz)));
            if (sgn == 0) {continue;}

            //------------------------------------------------------------------
            // History Update
            //------------------------------------------------------------------
            if (p < i-1)
            {
                // x[i] = x[i] * rhos[i-1]
                OK(SLIP_mpz_mul( SLIP_2D(x, i, k, mpz), SLIP_2D(x, i, k, mpz),
                                 rhos->x.mpz[i-1]));
                // x[i] = x[i] / rhos[p]
                if (p > -1)
                {
                    OK(SLIP_mpz_divexact( SLIP_2D(x, i, k, mpz), SLIP_2D(x, i, k, mpz),
                                          rhos->x.mpz[p]));
                }
            }

            //------------------------------------------------------------------
            // IPGE updates
            //------------------------------------------------------------------
            // Access the Lmi
            for (m = L->p[i]; m < L->p[i+1]; m++)
            {
                // Location of Lmi
                mnew = L->i[m];
                // skip if Lx[m] is zero
                OK(SLIP_mpz_sgn(&sgn, L->x.mpz[m]));
                if (sgn == 0) {continue;}
                // m > i
                if (mnew > i)
                {
                    p = h[mnew][k];
                    // x[mnew] is zero
                    OK(SLIP_mpz_sgn(&sgn, SLIP_2D(x, mnew, k, mpz)));
                    if (sgn == 0)
                    {
                        // x[m] = x[m] - lmi xi
                        OK(SLIP_mpz_submul( SLIP_2D(x, mnew, k, mpz), L->x.mpz[m],
                                            SLIP_2D(x, i, k, mpz)));
                        // x[m] = x[m] / rhos[i-1]
                        if (i > 0)
                        {
                            OK(SLIP_mpz_divexact( SLIP_2D(x, mnew, k, mpz), SLIP_2D(x, mnew, k, mpz), rhos->x.mpz[i-1]));
                        }
                    }
                    else
                    {
                        // History update if necessary
                        if (p < i-1)
                        {
                            // x[m] = x[m] * rhos[i-1]
                            OK(SLIP_mpz_mul( SLIP_2D(x, mnew, k, mpz), SLIP_2D(x, mnew, k, mpz),
                                             rhos->x.mpz[i-1]));
                            // x[m] = x[m] / rhos[p]
                            if (p > -1)
                            {
                                OK(SLIP_mpz_divexact( SLIP_2D(x, mnew, k, mpz), SLIP_2D(x, mnew, k, mpz),
                                                      rhos->x.mpz[p]));
                            }
                        }
                        // x[m] = x[m] * rhos[i]
                        OK(SLIP_mpz_mul( SLIP_2D(x, mnew, k, mpz), SLIP_2D(x, mnew, k, mpz),
                                         rhos->x.mpz[i]));
                        // x[m] = x[m] - lmi xi
                        OK(SLIP_mpz_submul( SLIP_2D(x, mnew, k, mpz), L->x.mpz[m], SLIP_2D(x, i, k, mpz)));
                        // x[m] = x[m] / rhos[i-1]
                        if (i > 0)
                        {
                            OK(SLIP_mpz_divexact( SLIP_2D(x, mnew, k, mpz), SLIP_2D(x, mnew, k, mpz), 
                                                  rhos->x.mpz[i-1]));
                        }
                    }
                    h[mnew][k] = i;
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // Free h memory
    //--------------------------------------------------------------------------
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

