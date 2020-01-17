//------------------------------------------------------------------------------
// IP_Chol/IP_Chol_permute_A: Symmetric permutation of matrix A
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* Purpose: Permute the matrix A and return A2 = PAP */
SLIP_sparse* IP_Chol_permute_A
(
    SLIP_sparse* A,        // Initial input matrix
    int* pinv,          // Row permutation
    SLIP_LU_analysis* S         // Column permutation
)
{
    SLIP_info ok;
    int nz = 0, j, n = A->n;
    SLIP_sparse* A2 = SLIP_create_sparse();
    OK(IP_sparse_alloc(A2, n, n, A->nz));
    if (ok != SLIP_OK)
        return NULL;
    OK(SLIP_mpq_set_ui(A2->scale, 1, 1));
    if (ok != SLIP_OK)
        return NULL;
    for (int k = 0 ; k < n ; k++)
    {
        A2->p [k] = nz ;                       // column k of A2 is column q[k] of A 
        j = S->q [k];
        for (int t = A->p [j] ; t < A->p [j+1] ; t++)
        {
            OK(SLIP_mpz_set(A2->x[nz], A->x[t]));  // row i of A is row pinv[i] of C 
            A2->i [nz++] = pinv [A->i [t]];
        }
    }
    A2->p [n] = nz ;                       // finalize the last column of C 
    A2->nz = A2->nzmax;
    return A2;
}