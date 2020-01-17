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
    mpq_t** x,              // rational solution to the system
    mpz_t** b,              // right hand side vector
    mpz_t* rhos,            // sequence of pivots
    SLIP_sparse* L,            // lower triangular matrix
    int* pinv,              // row permutation
    SLIP_options* option,// command options
    int numRHS              // number of RHS vectors
)
{
    SLIP_info ok;
    if (!x || !b || !rhos || !pinv) return SLIP_INCORRECT_INPUT;
    int i, k, n = L->n;
    // Permuted b
    mpz_t** b2 = SLIP_create_mpz_mat(n, numRHS);
    if (!b2) return SLIP_OUT_OF_MEMORY;
    // Set workspace b
    for (k = 0; k < numRHS; k++)
    {
        for (i = 0; i < n; i++)   
        {
            OK(SLIP_mpz_set(b2[pinv[i]][k], b[i][k]));
        }
    }
    // L*b2 = b2
    OK(IP_forward_sub(L, b2, rhos, numRHS));
    // b2 = b2 * det 
    OK(IP_array_mul(b2, rhos[n-1], n, numRHS));
    // U b2 = b2
    OK(IP_Chol_ltsolve(L, b2, numRHS));
    // x = b2/det 
    OK(IP_array_div(x, b2, rhos[n-1], n, numRHS));
    // Free memory
    FREE_WORKSPACE;
    return SLIP_OK;
}