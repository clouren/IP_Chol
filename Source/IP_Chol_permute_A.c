//------------------------------------------------------------------------------
// IP_Chol/IP_Chol_permute_A: Symmetric permutation of matrix A
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* Purpose: Permute the matrix A and return A2 = PAP */
SLIP_info IP_Chol_permute_A
(
    SLIP_matrix **A2_handle, // Output permuted matrix
    SLIP_matrix* A,          // Initial input matrix
    int64_t* pinv,           // Inverse row permutation
    SLIP_LU_analysis* S      // Column permutation
)
{
    SLIP_info ok;
    // Check inputs
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);
    if (!A2_handle || !pinv || !S)
        return SLIP_INCORRECT_INPUT;
    
    // Allocate memory for A2 which is a permuted copy of A
    int64_t nz = 0, j, n = A->n;
    SLIP_matrix* A2 = NULL;
    SLIP_matrix_allocate(&A2, SLIP_CSC, SLIP_MPZ, n, n, A->p[A->n], false, true, NULL);

    // Set A2 scale
    OK(SLIP_mpq_set(A2->scale, A->scale));
    // Populate the entries in A2
    for (int64_t k = 0 ; k < n ; k++)
    {
        A2->p [k] = nz ;                       // column k of A2 is column q[k] of A 
        j = S->q [k];
        for (int64_t t = A->p [j] ; t < A->p [j+1] ; t++)
        {
            OK(SLIP_mpz_set(A2->x.mpz[nz], A->x.mpz[t]));  // row i of A is row pinv[i] of A2 
            A2->i [nz++] = pinv [A->i [t]];
        }
    }
    A2->p [n] = nz ;                       // finalize the last column of A2 
    (*A2_handle) = A2;                     // Set A2_handle
    return SLIP_OK;
}
