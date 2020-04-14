//------------------------------------------------------------------------------
// IP_Chol/IP_Pre_Left_Factor: Symbolic left-looking Chol
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------


#include "../Include/IP-Chol.h"

/* Purpose: This function performs a symbolic left-looking factorization
 */
int IP_Pre_Left_Factor         // performs the Up looking Cholesky factorization
(
    SLIP_matrix* A,
    SLIP_matrix* L,              // partial L matrix
    int* xi,                  // nonzero pattern vector
    int* parent,              // Elimination tree
    Sym_chol * S,           // stores guess on nnz and column permutation
    int* c                   // Column pointers
)
{
    SLIP_info ok;
    // Input check/
    if (!L || !xi || !parent)
        return SLIP_INCORRECT_INPUT;
    
    int top, k, i, j, jnew, n = A->n;
    //--------------------------------------------------------------------------
    // Declare memory for L 
    //--------------------------------------------------------------------------
       
    // Allocate L  
    SLIP_matrix_allocate(&L, SLIP_CSC, SLIP_MPZ, n, n, S->lnz, false, false, option);
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
        
    L->i[0] = 0;
    c[0]++;

    //--------------------------------------------------------------------------
    // Iterations 1:n-1 (2:n in standard)
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
    {
        top = IP_Chol_ereach(A, k, parent, xi, c);  // Obtain nonzero pattern in xi[top..n]
     
        //----------------------------------------------------------------------
        // Iterate accross the nonzeros in x
        //----------------------------------------------------------------------
        int p = 0;
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
    L->nz = S->lnz;
    // Finalize L->p
    L->p[n] = L->nz;        
    return SLIP_OK;
}
