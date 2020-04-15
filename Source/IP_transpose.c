//------------------------------------------------------------------------------
// IP_Chol/IP_transpose: Transpose a matrix
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE  \
    SLIP_FREE(w);       \

#include "../Include/IP-Chol.h"
    
/* Purpose: This function sets C = A' */
SLIP_info IP_transpose
(
    SLIP_matrix **C_handle,     // C = A'
    SLIP_matrix *A              // Matrix to be transposed
)
{
    SLIP_info ok;
    int64_t* w = NULL;
    SLIP_matrix* C = NULL;
    SLIP_matrix_allocate(&C, SLIP_CSC, SLIP_MPZ, A->n, A->m, A->nz, false, true, NULL);
    int64_t p, q, j, n, m;
    m = A->m ; n = A->n ; 
    w = (int64_t*) SLIP_malloc(m* sizeof(int64_t));
    for (p = 0; p < m; p++) w[p] = 0;
    for (p = 0 ; p < A->p [n] ; p++) w [A->i [p]]++ ;       /* row counts */
    IP_cumsum_chol (C->p, w, m) ;                               /* row point64_t*ers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = A->p [j] ; p < A->p [j+1] ; p++)
        {
            q = w [A->i [p]]++;
            C->i [q] = j ;                 /* place A(i,j) as entry C(j,i) */
            OK(SLIP_mpz_set(C->x.mpz[q], A->x.mpz[p]));
        }
    }
    C->p[m] = A->p[n];
    (*C_handle) = C;
    FREE_WORKSPACE;
    return SLIP_OK;
}
