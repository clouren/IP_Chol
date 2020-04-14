//------------------------------------------------------------------------------
// IP_Chol/IP_Solve: Solve the linear system after factorization
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE                      \
    SLIP_delete_mpz_mat(&b2, n, numRHS);    \

#include "../Include/IP-Chol.h"
    
/* Purpose: This function solves the linear system LD^(-1)L' x = b.*/
int IP_Solve               //solves the linear system LD^(-1)L' x = b
(
    SLIP_matrix* x,              // rational solution to the system
    SLIP_matrix* b,              // right hand side vector
    SLIP_matrix* rhos,            // sequence of pivots
    SLIP_matrix* L,            // lower triangular matrix
    int* pinv,              // row permutation
    SLIP_options* option,// command options
    int numRHS              // number of RHS vectors
)
{
    SLIP_info ok;
    if (!x || !b || !rhos || !pinv) return SLIP_INCORRECT_INPUT;
    int i, k, n = L->n;
    // Permuted b
    SLIP_matrix *b2 = NULL;
    SLIP_matrix_copy(&b2, SLIP_DENSE SLIP_MPZ, b, option);
    for (k = 0; k < numRHS; k++)
    {
        for (i = 0; i < n; i++)   
        {
            OK (SLIP_mpz_set( SLIP_2D(b2, pinv[i], k, mpz),
                              SLIP_2D(b, i, k, mpz)));
        }
    }
    // L*b2 = b2
    OK(IP_forward_sub(L, b2, rhos, numRHS));
    // b2 = b2 * det 
    // Fix me
    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t k = 0; k < numRHS; k++)
        {
            // x[i][k] = x[i][k] * det
            OK(SLIP_mpz_mul(SLIP_2D(b2, i, k, mpz),
                            SLIP_2D(b2, i, k, mpz),
                            rhos->x.mpz[n-1]));
        }
    }
    
    
    // U b2 = b2
    OK(IP_Chol_ltsolve(L, b2, numRHS));
    // x = b2/det 
    mpq_t det2;
    mpq_init(det2);
    mpq_set_z(det2, rhos[n-1]);
    
    SLIP_matrix *x2 = NULL;
    SLIP_CHECK (SLIP_matrix_allocate(&x2, SLIP_DENSE, SLIP_MPQ, x->m, x->n,
        0, false, true, option)) ;
    
    //fixme
    for (int32_t i = 0; i < n; i++)
    {
        for (int32_t k = 0; k < numRHS; k++)
        {
            // Set x2[i] = x[i]
            OK(SLIP_mpq_set_num( SLIP_2D(x2, i, k, mpq), SLIP_2D(x, i, k, mpz));
            // x2[i] = x2[i] / det2
            OK(SLIP_mpq_div( SLIP_2D(x2, i, k, mpq), SLIP_2D(x2, i, k, mpq), det2));
        }
    }
    
    
    // Free memory
    FREE_WORKSPACE;
    return SLIP_OK;
}
