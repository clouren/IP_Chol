//------------------------------------------------------------------------------
// IP_Chol/IP_Solve: Solve the lint64_t*ear system after factorization
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE                      \
    SLIP_matrix_free(&b2, option);          \

#include "../Include/IP-Chol.h"
    
/* Purpose: This function solves the lint64_t*ear system LD^(-1)L' x = b.*/
SLIP_info IP_Solve               //solves the lint64_t*ear system LD^(-1)L' x = b
(
    // Output
    SLIP_matrix** x_handle,     // rational solution to the system
    // Input
    SLIP_matrix* A,             // Input matrix (permuted)
    SLIP_matrix* A_orig,        // Nonpermuted input matrix
    SLIP_matrix* b,             // right hand side vector
    SLIP_matrix* rhos,          // sequence of pivots
    SLIP_matrix* L,             // lower triangular matrix
    int64_t* pinv,                  // row permutation
    SLIP_LU_analysis *S,        // Column permutation
    SLIP_options* option        // command options
)
{
    SLIP_info ok;
    if (!x_handle || !b || !rhos || !pinv) return SLIP_INCORRECT_INPUT;
    
    int64_t i, j, k, n = L->n, nz;
    
    mpq_t scale ;
    SLIP_MPQ_SET_NULL (scale) ;

    SLIP_matrix *x = NULL;   // fint64_t*al solution
    SLIP_matrix *x2 = NULL;  // unpermuted solution
    SLIP_matrix *b2 = NULL;  // permuted b
    
    // Permuted b
    SLIP_matrix_allocate(&b2, SLIP_DENSE, SLIP_MPZ, b->m, b->n, b->m*b->n, false, true, option);
    for (i = 0; i < b->m; i++)
    {
        for (j = 0; j < b->n; j++)   
        {
            OK (SLIP_mpz_set( SLIP_2D(b2, pinv[i], j, mpz),
                              SLIP_2D(b, i, j, mpz)));
        }
    }
    
    // L*b2 = b2
    
    OK(IP_forward_sub(L, b2, rhos));
    
    // b2 = b2 * det 
    nz = SLIP_matrix_nnz(b, NULL);
    for (i = 0; i < nz; i++)
    {
        OK ( SLIP_mpz_mul( b2->x.mpz[i], b2->x.mpz[i], rhos->x.mpz[L->n-1]));
    }
    
    // U b2 = b2
    OK(IP_Chol_ltsolve(L, b2));
    
    // x = b2/det 
    mpq_t det2;
    mpq_init(det2);
    mpq_set_num(det2, rhos->x.mpz[L->n-1]);
    
    SLIP_CHECK (SLIP_matrix_allocate(&x, SLIP_DENSE, SLIP_MPQ, b2->m, b->n,
        0, false, true, option)) ;
        
    nz = SLIP_matrix_nnz(b, NULL);
    
    for (i = 0; i < nz; i++)
    {
        OK ( SLIP_mpq_set_num( x->x.mpq[i], b2->x.mpz[i]));
        OK ( SLIP_mpq_div( x->x.mpq[i], x->x.mpq[i], det2));
    }
    
    // Permute x
    //SLIP_matrix *x2 = NULL;
    SLIP_CHECK (SLIP_matrix_allocate(&x2, SLIP_DENSE, SLIP_MPQ, x->m, x->n,
        0, false, true, option)) ;
    
    for (int32_t i = 0; i < x->m; i++)
    {
        for (int32_t j = 0; j < x->n; j++)
        {
            SLIP_CHECK(SLIP_mpq_set( SLIP_2D(x2, S->q[i], j, mpq), SLIP_2D(x, i, j, mpq)));
        }
    }

    // Return x = x2
    for (int32_t i = 0; i < x->m; i++)
    {
        for (int32_t j = 0; j < x->n; j++)
        {
            SLIP_CHECK(SLIP_mpq_set( SLIP_2D(x, i, j, mpq), SLIP_2D(x2, i, j, mpq)));
        }
    }

    // Check solution
    bool check = option->check;
    if (check)
    {
        SLIP_CHECK (IP_check_solution (A_orig, x, b, option)) ;
    }
    
    
    //--------------------------------------------------------------------------
    // Scale the solution if necessary.
    //--------------------------------------------------------------------------

    SLIP_CHECK(SLIP_mpq_init(scale));

    
    
    // set the scalint64_t*g factor scale = A->scale / b->scale
    SLIP_CHECK( SLIP_mpq_div(scale, A->scale, b->scale));

    // Determine if the scalint64_t*g factor is 1
    int r;
    SLIP_CHECK(SLIP_mpq_cmp_ui(&r, scale, 1, 1));
    nz = x->m * x->n;
    if (r != 0 )
    {
        for (i = 0; i < nz; i++)
        {
            SLIP_CHECK(SLIP_mpq_mul(x->x.mpq[i], x->x.mpq[i], scale));
        }
    }
    
    
    // Free memory
    (*x_handle) = x;
    FREE_WORKSPACE;
    return SLIP_OK;
}
