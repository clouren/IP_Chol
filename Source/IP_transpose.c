//------------------------------------------------------------------------------
// REF_Chol/IP_transpose: Transpose a matrix
//------------------------------------------------------------------------------

// REF Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// Texas A&M University.  All Rights Reserved.  See REF_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE  \
    SLIP_FREE(w);       \

#include "../Include/REF-Chol.h"
    
/* Purpose: This function sets C = A' 
 * C_handle is NULL on input. On output, C_handle contains a pointer to A'
 */
IP_Chol_info IP_transpose
(
    SLIP_matrix **C_handle,     // C = A'
    SLIP_matrix *A              // Matrix to be transposed
)
{
    IP_Chol_info ok;
    // Check input
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);
    if (!C_handle) 
        return SLIP_INCORRECT_INPUT;
    
    // Declare workspace and C
    int64_t* w = NULL;
    SLIP_matrix* C = NULL;
    OK(SLIP_matrix_allocate(&C, SLIP_CSC, SLIP_MPZ, A->n, A->m, A->p[A->n], false, true, NULL));
    int64_t p, q, j, n, m;
    m = A->m ; n = A->n ; 
    w = (int64_t*) SLIP_malloc(m* sizeof(int64_t));
    if (!w)
    {
        SLIP_matrix_free(&C, NULL);
        return SLIP_OUT_OF_MEMORY;
    }
    for (p = 0; p < m; p++) w[p] = 0;
    for (p = 0 ; p < A->p [n] ; p++) w [A->i [p]]++ ;       // row counts 
    IP_cumsum_chol (C->p, w, m) ;                           // row pointers
    for (j = 0 ; j < n ; j++)
    {
        for (p = A->p [j] ; p < A->p [j+1] ; p++)
        {
            q = w [A->i [p]]++;
            C->i [q] = j ;                 // place A(i,j) as entry C(j,i) 
            OK(SLIP_mpz_set(C->x.mpz[q], A->x.mpz[p]));
        }
    }
    C->p[m] = A->p[n];
    (*C_handle) = C;
    FREE_WORKSPACE;
    return SLIP_OK;
}
