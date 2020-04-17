//------------------------------------------------------------------------------
// IP_Chol/IP_Solve: Solve the lint64_t*ear system after factorization
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE                      \
    SLIP_matrix_free(&b2, option);          \
    SLIP_matrix_free(&x, option);           \
    SLIP_MPQ_CLEAR(scale);                  \
    SLIP_MPQ_CLEAR(det2);                   \

#include "../Include/IP-Chol.h"
    
/* Purpose: This function solves the linear system LD^(-1)L' x = b.*
 *
 * Input arguments:
 * 
 * x_handle:        A handle to the solution matrix. On input this is NULL,
 *                  on output x_handle contains a pointer to the solution vector(s)
 * 
 * A:               Permuted version of the input matrix
 * 
 * A_orig:          Nonpermuted input matrix
 * 
 * b:               Right hand side vector(s)
 * 
 * rhos:            Sequence of pivots encountered in factorization
 * 
 * L:               Lower triangular matrix
 * 
 * pinv:            Inverse row permutation
 * 
 * S:               Column permutation struct
 * 
 * option:          Command options
 * 
 */
SLIP_info IP_Solve              // solves the linear system LDL' x = b
(
    // Output
    SLIP_matrix** x_handle,     // rational solution to the system
    // Input
    SLIP_matrix* A,             // Input matrix (permuted)
    SLIP_matrix* A_orig,        // Nonpermuted input matrix
    SLIP_matrix* b,             // right hand side vector
    SLIP_matrix* rhos,          // sequence of pivots
    SLIP_matrix* L,             // lower triangular matrix
    int64_t* pinv,              // row permutation
    SLIP_LU_analysis *S,        // Column permutation
    SLIP_options* option        // command options
)
{
    SLIP_info ok;
    // Check the inputs
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);
    SLIP_REQUIRE(A_orig, SLIP_CSC, SLIP_MPZ);
    SLIP_REQUIRE(b, SLIP_DENSE, SLIP_MPZ);
    SLIP_REQUIRE(rhos, SLIP_DENSE, SLIP_MPZ);
    SLIP_REQUIRE(L, SLIP_CSC, SLIP_MPZ);
    
    if (!x_handle || !pinv || !S) return SLIP_INCORRECT_INPUT;
    
    int64_t i, j, k, n = L->n, nz;
    
    mpq_t scale, det2 ;
    SLIP_MPQ_SET_NULL (scale) ;
    SLIP_MPQ_SET_NULL (det2) ;

    SLIP_matrix *x = NULL;   // unpermuted solution
    SLIP_matrix *x2 = NULL;  // permuted final solution
    SLIP_matrix *b2 = NULL;  // permuted b
    
    // Permute b and place it in b2
    SLIP_matrix_allocate(&b2, SLIP_DENSE, SLIP_MPZ, b->m, b->n, b->m*b->n, false, true, option);
    for (i = 0; i < b->m; i++)
    {
        for (j = 0; j < b->n; j++)   
        {
            OK (SLIP_mpz_set( SLIP_2D(b2, pinv[i], j, mpz),
                              SLIP_2D(b, i, j, mpz)));
        }
    }
    
    // b2 = L \ b2    
    OK(IP_forward_sub(L, b2, rhos));
    
    // b2 = b2 * det 
    nz = SLIP_matrix_nnz(b, NULL);
    for (i = 0; i < nz; i++)
    {
        OK ( SLIP_mpz_mul( b2->x.mpz[i], b2->x.mpz[i], rhos->x.mpz[L->n-1]));
    }
    
    SLIP_CHECK(SLIP_mpq_init(det2));
    // b2 = L' \ b2
    OK(IP_Chol_ltsolve(L, b2));
    
    // x = b2/det 
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
    SLIP_CHECK (SLIP_matrix_allocate(&x2, SLIP_DENSE, SLIP_MPQ, x->m, x->n,
        0, false, true, option)) ;
    
    for (int32_t i = 0; i < x->m; i++)
    {
        for (int32_t j = 0; j < x->n; j++)
        {
            SLIP_CHECK(SLIP_mpq_set( SLIP_2D(x2, S->q[i], j, mpq), SLIP_2D(x, i, j, mpq)));
        }
    }

    // Check solution
    bool check = option->check;
    if (check)
    {
        SLIP_CHECK (IP_check_solution (A_orig, x2, b, option)) ;
    }
    
    
    //--------------------------------------------------------------------------
    // Scale the solution if necessary.
    //--------------------------------------------------------------------------

    SLIP_CHECK(SLIP_mpq_init(scale));

    // set the scaling factor scale = A->scale / b->scale
    SLIP_CHECK( SLIP_mpq_div(scale, A->scale, b->scale));

    // Determine if the scaling factor is 1
    int r;
    SLIP_CHECK(SLIP_mpq_cmp_ui(&r, scale, 1, 1));
    nz = x->m * x->n;
    if (r != 0 )
    {
        for (i = 0; i < nz; i++)
        {
            SLIP_CHECK(SLIP_mpq_mul(x2->x.mpq[i], x2->x.mpq[i], scale));
        }
    }
    
    
    // Free memory
    (*x_handle) = x2;
    FREE_WORKSPACE;
    return SLIP_OK;
}
