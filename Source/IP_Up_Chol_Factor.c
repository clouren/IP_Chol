//------------------------------------------------------------------------------
// IP_Chol/IP_Up_Chol_Factor: Up-lookint64_t*g Cholesky factorization
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE              \
    SLIP_matrix_free(&x, NULL);     \
    SLIP_FREE(c);                   \
    SLIP_FREE(xi);                  \
    SLIP_FREE(h);                   \
    SLIP_FREE(post);                \


#include "../Include/IP-Chol.h"
    
/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
SLIP_info IP_Up_Chol_Factor           // performs the Up lookint64_t*g Cholesky factorization
(
    SLIP_matrix* A,             // matrix to be factored
    SLIP_matrix** L_handle,     // lower triangular matrix
    Sym_chol * S,               // stores guess on nnz and column permutation
    SLIP_matrix ** rhos_handle, // sequence of pivots
    SLIP_options* option        // command options
)
{
    SLIP_info ok;
    // Input check
    if (!A || !L_handle || !S || !rhos_handle || !option )
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    int64_t anz = SLIP_matrix_nnz (A, option) ;

    if (anz < 0)
    {
        return SLIP_INCORRECT_INPUT;
    }

    (*L_handle) = NULL ;
    (*rhos_handle) = NULL ;
    
    
    //--------------------------------------------------------------------------
    // Declare and int64_t*itialize workspace 
    //--------------------------------------------------------------------------
    
    SLIP_matrix *L = NULL ;
    SLIP_matrix *rhos = NULL ;
    int64_t *xi = NULL ;
    int64_t *h = NULL ;
    SLIP_matrix *x = NULL ;
    
    //// Begint64_t* timint64_t*g factorization
    int64_t n = A->n, top, i, j, col, loc, lnz = 0, unz = 0, pivot, jnew, k;
    size_t size;

    int64_t* post = NULL;
    int64_t* c = NULL;

    
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

    // Obtaint64_t* the elimint64_t*ation tree of A
    S->parent = IP_Chol_etree(A);              // Obtaint64_t* the elim tree
    post = IP_Chol_post(S->parent, n);    // Postorder the tree
    
    // Get the column counts of A
    c = IP_Chol_counts(A, S->parent, post, 0);
    
    S->cp = (int64_t*) SLIP_malloc( (n+1)*sizeof(int64_t*));
    S->lnz = IP_cumsum_chol(S->cp, c, n);    // Get column point64_t*ers for L
   
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
    // Declare memory for L 
    //--------------------------------------------------------------------------
    
    
    // Allocate L without int64_t*itializint64_t*g each entry.
    // L is allocated to have nnz(L) which is estimated by the symbolic
    // analysis. However, unlike traditional matrix allocation, the second
    // boolean parameter here is set to false, so the int64_t*dividual values of
    // L are not allocated. Instead, a more efficient method to
    // allocate these values is done int64_t* the factorization to reduce
    // memory usage.
    SLIP_CHECK (SLIP_matrix_allocate(&L, SLIP_CSC, SLIP_MPZ, n, n, S->lnz,
        false, false, option));
    
    // Set the column  point64_t*ers of L
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
    

    //--------------------------------------------------------------------------
    // Iterations 0:n-1 (1:n int64_t* standard)
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        // LDx = A(:,k)
        SLIP_CHECK ( IP_Up_Chol_triangular_solve(&top, L, A, k, xi, S->parent, c, rhos, h, x));
        
        // If x[k] is nonzero that is the pivot. if x[k] == 0 then matrix is sint64_t*gular.
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
        int64_t p = 0;
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
        OK(SLIP_mpz_set(L->x.mpz[p], x->x.mpz[k]));
    }
    // Finalize L->p
    L->p[n] = S->lnz;        

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    (*L_handle) = L;
    (*rhos_handle) = rhos;
    FREE_WORKSPACE;
    return SLIP_OK;
}
