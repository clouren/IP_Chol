//------------------------------------------------------------------------------
// IP_Chol/IP_check_solution: check solution to Ax=b
//------------------------------------------------------------------------------
//
//
// IP Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// and Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// IP_Chol/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Check the solution of the linear system by performing a rational
 * arithmetic A*x and checking if (A*x) == b. This function should only be 
 * used for debugging purposes.
 */

#define SLIP_FREE_ALL                       \
    SLIP_MPQ_CLEAR(temp);                   \
    SLIP_matrix_free(&b2, NULL);

#include "../Include/IP-Chol.h"

IP_Chol_info IP_check_solution
(
    const SLIP_matrix *A,         // Input matrix
    const SLIP_matrix *x,         // Solution vectors
    const SLIP_matrix *b,         // Right hand side vectors
    const SLIP_options* option    // Command options
)
{

    //--------------------------------------------------------------------------
    // check int64_t*puts
    //--------------------------------------------------------------------------

    IP_Chol_info ok, info = IP_Chol_OK ;
    SLIP_REQUIRE (A, SLIP_CSC,   SLIP_MPZ) ;
    SLIP_REQUIRE (x, SLIP_DENSE, SLIP_MPQ) ;
    SLIP_REQUIRE (b, SLIP_DENSE, SLIP_MPZ) ;

    //--------------------------------------------------------------------------
    // Declare vars
    //--------------------------------------------------------------------------

    int64_t p, j, i ;
    SLIP_matrix *b2 = NULL;   // b2 stores the solution of A*x
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);

    OK (SLIP_mpq_init(temp));
    OK (SLIP_matrix_allocate(&b2, SLIP_DENSE, SLIP_MPQ, b->m, b->n,
        b->nzmax, false, true, option));

    //--------------------------------------------------------------------------
    // perform SLIP_mpq_addmul in loops
    //--------------------------------------------------------------------------

    for (j = 0; j < b->n; j++)
    {
        for (i = 0; i < b->m; i++)
        {
            for (p = A->p[i]; p < A->p[i + 1]; p++)
            {
                // temp = A[p][i]
                OK(SLIP_mpq_set_z(temp, A->x.mpz[p]));

                // temp = temp*x[i]
                OK(SLIP_mpq_mul(temp, temp,
                                        SLIP_2D(x, i, j, mpq)));

                // b2[p] = b2[p]-temp
                OK(SLIP_mpq_add(SLIP_2D(b2, A->i[p], j, mpq),
                                        SLIP_2D(b2, A->i[p], j, mpq),temp));
            }
        }
    }

    //--------------------------------------------------------------------------
    // check if b==b2
    //--------------------------------------------------------------------------

    for (j = 0; j < b->n; j++)
    {
        for (i = 0; i < b->m; i++)
        {
            // z = b[i] (correct b)
            OK (SLIP_mpq_set_z(temp, SLIP_2D(b, i, j, mpz)));

            // set check false if b!=b2
            int r ;
            OK (SLIP_mpq_equal(&r, temp, SLIP_2D(b2, i, j, mpq)));
            if (r == 0)
            {
                info = SLIP_INCORRECT;
                j = b->n;
                break;
            }
        }
    }

    //--------------------------------------------------------------------------
    // Print info if necessary
    //--------------------------------------------------------------------------

    int64_t pr = option->print_level;
    if (pr >= 1)
    {
        if (info == IP_Chol_OK)
        {
            printf ("Solution is verified to be exact.\n") ;
        }
        else if (info == IP_Chol_INCORRECT)
        {
            // This can never happen.
            printf ("ERROR! Solution is wrong. This is a bug; please "
                    "contact the authors of IP_Chol.\n") ;
        }
    }

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    SLIP_FREE_ALL;
    return info;
}

