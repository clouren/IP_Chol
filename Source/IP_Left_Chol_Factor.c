//------------------------------------------------------------------------------
// IP_Chol/IP_Left_Chol_Factor: Compute the left-lookint64_t*g cholesky factorization of SPD matrix A
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

//TODO Combine int o a single IP_Chol_factor function which has a bool int64_t*put which decides left or up

#define FREE_WORKSPACE              \
    SLIP_matrix_free(&x,option);    \
    SLIP_FREE(xi);                  \
    SLIP_FREE(h);                   \
    SLIP_FREE(c);                   \
    SLIP_FREE(post);                \
    
#include "../Include/IP-Chol.h"

static inline int64_t compare3 (const void * a, const void * b)
{
    return ( *(int64_t*)a - *(int64_t*)b );
}


/* Purpose: This function performs the Left-lookint64_t*g IP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
SLIP_info IP_Left_Chol_Factor         // performs the SLIP LU factorization
(
    SLIP_matrix* A,           // matrix to be factored
    SLIP_matrix** L_handle,           // lower triangular matrix
    Sym_chol * S,           // stores guess on nnz and column permutation
    SLIP_matrix ** rhos_handle, // sequence of pivots
    SLIP_options* option// command options
)
{
    SLIP_info ok;
    
    // Input check
    if (!A || !L_handle || !S || !rhos_handle || !option )
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    int64_t  anz = SLIP_matrix_nnz (A, option) ;

    if (anz < 0)
    {
        return SLIP_INCORRECT_INPUT;
    }

    (*L_handle) = NULL ;
    (*rhos_handle) = NULL ;
    
    //--------------------------------------------------------------------------
    // Declare and int64_t*itialize workspace 
    //--------------------------------------------------------------------------
    
    SLIP_matrix *rhos = NULL ;
    int64_t  *xi = NULL ;
    int64_t  *h = NULL ;
    SLIP_matrix *x = NULL ;
    
    //// Begint64_t* timint64_t*g factorization
    int64_t  n = A->n, top, i, j, k=0, col, loc, lnz = 0, unz = 0, pivot, jnew;
    size_t size;

    int64_t* post = NULL;
    int64_t* c = NULL;
    
    
    n = A->n;
    
    // h is the history vector utilized for the sparse REF
    // triangular solve algorithm. h serves as a global
    // vector which is repeatedly passed int64_t*o the triangular
    // solve algorithm
    h = (int64_t*) SLIP_malloc(n* sizeof(int64_t));
    
    // xi is the global nonzero pattern vector. It stores
    // the pattern of nonzeros of the kth column of L and U
    // for the triangular solve.
    xi = (int64_t*) SLIP_malloc(2*n* sizeof(int64_t));
    if (!h || !xi) 
    {
        FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }  
    
    // int64_t*itialize workspace history array
    for (i = 0; i < n; i++)
    {
        h[i] = -1;
    }
    
    // Obtaint64_t* the elimination tree of A
    S->parent = IP_Chol_etree(A);              // Obtain the elim tree
    post = IP_Chol_post(S->parent, n);    // Postorder the tree
    
    // Get the column counts of A
    c = IP_Chol_counts(A, S->parent, post, 0);
    
    S->cp = (int64_t*) SLIP_malloc( (n+1)* sizeof(int64_t));
    
    if (!S->parent || !post || !c || !S->cp)
    {
        FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    
    S->lnz = IP_cumsum_chol(S->cp, c, n);    // Get column pointers for L

   //--------------------------------------------------------------------------
    // allocate and int64_t*itialize the workspace x
    //--------------------------------------------------------------------------

    // SLIP LU utilizes arbitrary sized int64_t*egers which can grow beyond the
    // default 64 bits allocated by GMP. If the int64_t*egers frequently grow, GMP
    // can get bogged down by performint64_t*g int64_t*ermediate reallocations. Instead,
    // we utilize a larger estimate on the workspace x vector so that computint64_t*g
    // the values int64_t* L and U do not require too many extra int64_t*emediate calls to
    // realloc.
    //
    // Note that the estimate presented here is not an upper bound nor a lower
    // bound.  It is still possible that more bits will be required which is
    // correctly handled int64_t*ernally.
    int64_t estimate = 64 * SLIP_MAX (2, ceil (log2 ((double) n))) ;

    // Create x, a global dense mpz_t matrix of dimension n*1. Unlike rhos, the
    // second boolean parameter is set to false to avoid int64_t*itializint64_t*g
    // each mpz entry of x with default size.  It is int64_t*ialized below.
    SLIP_CHECK (SLIP_matrix_allocate(&x, SLIP_DENSE, SLIP_MPZ, n, 1, n,
        false, /* do not int64_t*itialize the entries of x: */ false, option));
    
    // Create rhos, a global dense mpz_t matrix of dimension n*1. 
    SLIP_CHECK (SLIP_matrix_allocate(&rhos, SLIP_DENSE, SLIP_MPZ, n, 1, n,
        false, true, option));
    
    if (!x)
    {
        FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    
    // int64_t*itialize the entries of x
    for (i = 0; i < n; i++)
    {
        // Allocate memory for entries of x
        SLIP_CHECK(SLIP_mpz_init2(x->x.mpz[i], estimate));
    }

    //--------------------------------------------------------------------------
    // Declare memory for x, L 
    //--------------------------------------------------------------------------
    
    // Allocate L  
    SLIP_matrix * L = NULL;
    OK(IP_Pre_Left_Factor(A, &L, xi, S->parent, S, c));
        
    
    //--------------------------------------------------------------------------
    // Set C
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        L->p[k] = c[k] = L->p[k];
    }
    
  
  
  //--------------------------------------------------------------------------
    // Iterations 0:n-1 (2:n int64_t* standard)
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        // LDx = A(:,k)
        OK (IP_Left_Chol_triangular_solve(&top, L, A, k, xi, rhos, h, x, S->parent, c));

        if (mpz_sgn(x->x.mpz[k]) != 0)
        {
            pivot = k;
            OK(SLIP_mpz_set(rhos->x.mpz[k], x->x.mpz[k]));
        }
        else
        {
            FREE_WORKSPACE;
            return SLIP_SINGULAR;
        }
        //----------------------------------------------------------------------
        // Iterate accross the nonzeros int64_t* x
        //----------------------------------------------------------------------
        for (j = top; j < n; j++)
        {
            jnew = xi[j];
            if (jnew >= k)
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
        //printf("\nEnd of iteration %ld L is:\n", k);
        //SLIP_matrix_check(L,option);
    }
    L->p[n] = S->lnz;

    //--------------------------------------------------------------------------
    // Free memory
    
    (*L_handle) = L;
    (*rhos_handle) = rhos;
    FREE_WORKSPACE;
    return SLIP_OK;
}
