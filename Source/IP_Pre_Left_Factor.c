//------------------------------------------------------------------------------
// IP_Chol/IP_Pre_Left_Factor: Symbolic left-looking Chol
//------------------------------------------------------------------------------

// IP Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// and Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// IP_Chol/License for the license.

//------------------------------------------------------------------------------


#include "../Include/IP-Chol.h"

/* Purpose: This function performs a symbolic left-looking factorization
 * It allocates the memory for the L matrix and allocates the individual
 * entries in the matrix.
 */
IP_Chol_info IP_Pre_Left_Factor         // performs a symbolic Cholesky factorization
(
    SLIP_matrix* A,                 // Input matrix
    SLIP_matrix** L_handle,         // partial L matrix
    int64_t*  xi,                   // nonzero pattern vector
    int64_t*  parent,               // Elimination tree
    Sym_chol * S,                   // stores nnz and elimination tree
    int64_t * c                     // Column point64_t*ers
)
{
    // Input check/
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);
    if (!L_handle || !xi || !parent || !S || !c)
        return SLIP_INCORRECT_INPUT;
    
    int64_t  top, k, j, jnew, n = A->n;
    //--------------------------------------------------------------------------
    // Declare memory for L 
    //--------------------------------------------------------------------------
       
    // Allocate L  
    SLIP_matrix* L = NULL;
    SLIP_matrix_allocate(&L, SLIP_CSC, SLIP_MPZ, n, n, S->lnz, false, false, NULL);
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
        
    L->i[0] = 0;
    c[0]++;

    //--------------------------------------------------------------------------
    // Iterations 1:n-1 (2:n int64_t* standard)
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
    {
        top = IP_Chol_ereach(A, k, parent, xi, c);  // Obtain nonzero pattern in xi[top..n]
     
        //----------------------------------------------------------------------
        // Iterate accross the nonzeros in x
        //----------------------------------------------------------------------
        int64_t p = 0;
        for (j = top; j < n; j++)
        {
            jnew = xi[j];
            if (jnew == k) continue;
            p = c[jnew]++;
            // Place the i location of the L->nz nonzero
            L->i[p] = k;
        }
        p = c[k]++;
        L->i[p] = k;
    }
    // Finalize L->p
    L->p[n] = S->lnz;
    (*L_handle) = L;
    return SLIP_OK;
}
