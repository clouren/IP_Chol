//------------------------------------------------------------------------------
// IP_Chol/IP_Up_Chol_Factor: Up-looking Cholesky factorization
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE              \
    SLIP_delete_mpz_array(&x,n);    \
    SLIP_FREE(c);                   \
    SLIP_FREE(xi);                  \
    SLIP_FREE(h);                   \
    SLIP_FREE(post);                \
    SLIP_MPFR_CLEAR(temp);          \
    SLIP_MPZ_CLEAR(sigma);          \


#include "../Include/IP-Chol.h"
    
/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
int IP_Up_Chol_Factor         // performs the Up looking Cholesky factorization
(
    SLIP_matrix* A,           // matrix to be factored
    SLIP_matrix** L_handle,           // lower triangular matrix
    Sym_chol * S,           // stores guess on nnz and column permutation
    SLIP_matrix ** rhos_handle, // sequence of pivots
    int** pinv_handle,             // inverse row permutation
    SLIP_options* option// command options
)
{
    SLIP_info ok;
    // Input check
    if (!A || !L || !S || !rhos || !option)
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    int64_t anz = SLIP_matrix_nnz (A, option) ;

    if (!L_handle  || !rhos_handle || !pinv_handle || !S || anz < 0)
    {
        return SLIP_INCORRECT_INPUT;
    }

    (*L_handle) = NULL ;
    (*rhos_handle) = NULL ;
    (*pinv_handle) = NULL ;
    
    
    //--------------------------------------------------------------------------
    // Declare and initialize workspace 
    //--------------------------------------------------------------------------
    //// Begin timing factorization
    int32_t n = A->n, numRealloc = 0, k = 0, top, i, j, col, loc, lnz = 0, unz = 0, pivot,
        check, jnew;
    long size;

    int* h = NULL;
    int* xi = NULL;
    int32_t* post = NULL;
    int32_t* c = NULL;
    SLIP_matrix* x = NULL;

    
    // History vector
    h = (int32_t*) SLIP_malloc(n* sizeof(int32_t));
    // Nonzero pattern
    xi = (int32_t*) SLIP_malloc(2*n* sizeof(int32_t));

    if (!h || !xi)
    {
        FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    IP_reset_int_array(h,n);

    // Obtain the elimination tree of A
    S->parent = IP_Chol_etree(A);              // Obtain the elim tree
    post = IP_Chol_post(S->parent, n);    // Postorder the tree
    
    // Get the column counts of A
    c = IP_Chol_counts(A, S->parent, post, 0);
    
    S->cp = (int32_t*) SLIP_malloc( (n+1)*sizeof(int32_t));
    S->lnz = IP_cumsum_chol(S->cp, c, n);    // Get column pointers for L
   
     //--------------------------------------------------------------------------
    // allocate and initialize the workspace x
    //--------------------------------------------------------------------------

    // SLIP LU utilizes arbitrary sized integers which can grow beyond the
    // default 64 bits allocated by GMP. If the integers frequently grow, GMP
    // can get bogged down by performing intermediate reallocations. Instead,
    // we utilize a larger estimate on the workspace x vector so that computing
    // the values in L and U do not require too many extra intemediate calls to
    // realloc.
    //
    // Note that the estimate presented here is not an upper bound nor a lower
    // bound.  It is still possible that more bits will be required which is
    // correctly handled internally.
    int64_t estimate = 64 * SLIP_MAX (2, ceil (log2 ((double) n))) ;

    // Create x, a global dense mpz_t matrix of dimension n*1. Unlike rhos, the
    // second boolean parameter is set to false to avoid initializing
    // each mpz entry of x with default size.  It is intialized below.
    SLIP_CHECK (SLIP_matrix_allocate(&x, SLIP_DENSE, SLIP_MPZ, n, 1, n,
        false, /* do not initialize the entries of x: */ false, option));
    
    

    //--------------------------------------------------------------------------
    // Declare memory for x, L 
    //--------------------------------------------------------------------------
    // Initialize x
    if (!x)
    {
        FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    
    // Allocate L  
    SLIP_matrix * L = NULL;
    SLIP_matrix_allocate(&L, SLIP_CSC, SLIP_MPZ, n, n, S->lnz, false. false, option);
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
    

    //--------------------------------------------------------------------------
    // Iterations 0:n-1 (1:n in standard)
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        // LDx = A(:,k)
        top = IP_Up_Chol_triangular_solve(L, A, k, xi, S->parent, c, rhos, h, x);
        if (top < 0)
        {
            FREE_WORKSPACE;
            return top;
        }
        if (mpz_sgn(x[k]) != 0)
        {
            pivot = k;
            OK(SLIP_mpz_set(rhos[k], x[k]));
        }
        else
        {
            pivot = SLIP_SINGULAR;
            FREE_WORKSPACE;
            return pivot;
        }
        
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
            size = mpz_sizeinbase(x->x.mpz[jnew],2);
            // GMP manual: Allocated size should be size+2
            OK(SLIP_mpz_init2(L->x.mpz[p], size+2));
            // Place the x value of the L->nz nonzero
            OK(SLIP_mpz_set(L->x.mpz[p],x->x.mpz[jnew]));
        }
        p = c[k]++;
        L->i[p] = k;
        size = mpz_sizeinbase(x->x.mpz[k], 2);
        OK(SLIP_mpz_init2(L->x.mpz[p], size+2));
        OK(SLIP_mpz_set(L->x.mpz[p], x->x.mpx[k]));
    }
    L->nz = S->lnz;
    // Finalize L->p
    L->p[n] = L->nz;        

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
    return SLIP_OK;
}
