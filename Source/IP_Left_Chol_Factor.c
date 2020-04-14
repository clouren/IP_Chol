//------------------------------------------------------------------------------
// IP_Chol/IP_Left_Chol_Factor: Compute the left-looking cholesky factorization of SPD matrix A
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE              \
    SLIP_delete_mpz_array(&x,n);    \
    SLIP_FREE(xi);                  \
    SLIP_FREE(h);                   \
    SLIP_FREE(c);                   \
    SLIP_FREE(post);                \
    SLIP_MPFR_CLEAR(temp);          \
    SLIP_MPZ_CLEAR(sigma);          \
    
#include "../Include/IP-Chol.h"

static inline int32_t compare3 (const void * a, const void * b)
{
    return ( *(int32_t*)a - *(int32_t*)b );
}


/* Purpose: This function performs the Left-looking IP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
int IP_Left_Chol_Factor         // performs the SLIP LU factorization
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
    if (!A || !L || !S || !rhos || !pinv || !option)
        return SLIP_INCORRECT_INPUT;
    //--------------------------------------------------------------------------
    // Declare and initialize workspace 
    //--------------------------------------------------------------------------
    int32_t n, numRealloc = 0, k = 0, top, i, j, col, loc, lnz = 0, unz = 0, pivot,
        check, jnew;
    long size;

    
    int64_t anz = SLIP_matrix_nnz (A, option) ;

    if (!L_handle  || !rhos_handle || !pinv_handle || !S || anz < 0)
    {
        return SLIP_INCORRECT_INPUT;
    }

    (*L_handle) = NULL ;
    (*rhos_handle) = NULL ;
    (*pinv_handle) = NULL ;
    
    
    
    int32_t* h = NULL;
    int32_t* xi = NULL;
    int32_t* c = NULL;
    int32_t* post = NULL;
    SLIP_matrix* x = NULL;
    
    
    n = A->n;
    // History vector
    h = (int32_t*) SLIP_malloc(n*sizeof(int32_t));
    // Nonzero pattern
    xi = (int32_t*) SLIP_malloc(2*n*sizeof(int32_t));
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
    
    S->cp = (int32_t*) SLIP_malloc( (n+1)* sizeof(int32_t));
    
    if (!S->parent || !post || !c || !S->cp)
    {
        FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    
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
    for (i = 0; i < n; i++) pinv[i] = i;
    
    // Allocate L  
    SLIP_matrix * L = NULL;
    OK(IP_Pre_Left_Factor(A, &L, xi, S->parent, S, c));
    
    //--------------------------------------------------------------------------
    // Iteration 0, must select pivot
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        c[k] = L->p[k];
    }
  
    //--------------------------------------------------------------------------
    // Iterations 0:n-1 (2:n in standard)
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        // LDx = A(:,k)
        top = IP_Left_Chol_triangular_solve(L, A, k, xi, rhos, pinv, h, x, S->parent, c);
        if (top < 0)
        {
            FREE_WORKSPACE;
            return top;
        }
        if (mpz_sgn(x->x.mpz[k]) != 0)
        {
            pivot = k;
            OK(SLIP_mpz_set(rhos[k], x[k]));
        }
        else
            pivot = SLIP_SINGULAR;
            
        // Error
        if (pivot < SLIP_OK)
        {
            FREE_WORKSPACE;
            break;
        }
        
        //----------------------------------------------------------------------
        // Iterate accross the nonzeros in x
        //----------------------------------------------------------------------
        for (j = top; j < n; j++)
        {
            jnew = xi[j];
            if (jnew >= k && ok == SLIP_OK)
            {
                // Place the i location of the L->nz nonzero
                size = mpz_sizeinbase(x->x.mpz[jnew],2);
                // GMP manual: Allocated size should be size+2
                OK(SLIP_mpz_init2(L->x.mpz[lnz], size+2));
                // Place the x value of the L->nz nonzero
                OK(SLIP_mpz_set(L->x.mpz[lnz],x->x.mpz[jnew]));
                // Increment L->nz
                lnz += 1;
            }
        }
    }
    L->nz = S->lnz;

    //--------------------------------------------------------------------------
    // Free memory
    FREE_WORKSPACE;
    return ok;
}
