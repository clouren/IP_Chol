//------------------------------------------------------------------------------
// REF_Chol/IP_Chol_Factor: Integer preserving Cholesky factorization
//------------------------------------------------------------------------------

// REF Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// Texas A&M University.  All Rights Reserved.  See REF_Chol/License for the license.

//------------------------------------------------------------------------------

#define FREE_WORKSPACE              \
    SLIP_matrix_free(&x, NULL);     \
    SLIP_FREE(c);                   \
    SLIP_FREE(xi);                  \
    SLIP_FREE(h);                   \
    SLIP_FREE(post);                \


#include "../Include/REF-Chol.h"
    
/* Purpose: This function performs the integer preserving Cholesky factorization.
 * It allows either the left-looking or up-looking integer-preserving Cholesky factorization.
 * In order to compute the L matrix, it performs n iterations of a sparse REF symmetric
 * triangular solve function. The overall factorization is PAP' = LDL'
 * 
 * Input arguments:
 * 
 * A:           The user's permuted input matrix
 * 
 * L_handle:    A handle to the L matrix. Null on input. On output, contains a pointer to the 
 *              L matrix
 * 
 * S:           Symbolic analysis struct for Cholesky factorization. NULL on input. On output,
 *              contains the elimination tree and number of nonzeros in L.
 * 
 * rhos_handle: A handle to the sequence of pivots. NULL on input. On output, contains a pointer
 *              to the pivots matrix.
 * 
 * left:        A boolean parameter which tells the function whether it is performing a left-looking
 *              or up-looking factorization. If this bool is true, a left-looking factorization
 *              is done, otherwise the up-looking factorization is done.
 * 
 * option:      Command options
 * 
 */
IP_Chol_info IP_Chol_Factor        // performs an integer-preserving Cholesky factorization
(
    SLIP_matrix* A,             // matrix to be factored
    SLIP_matrix** L_handle,     // lower triangular matrix
    Sym_chol * S,               // stores guess on nnz and column permutation
    SLIP_matrix ** rhos_handle, // sequence of pivots
    bool left,                  // Set true if performing a left-looking factorization
    SLIP_options* option        // command options
)
{
    IP_Chol_info ok;
    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);
    if (!L_handle || !S || !rhos_handle || !option )
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
    // Declare and initialize workspace 
    //--------------------------------------------------------------------------
    
    SLIP_matrix *L = NULL ;
    SLIP_matrix *rhos = NULL ;
    int64_t *xi = NULL ;
    int64_t *h = NULL ;
    SLIP_matrix *x = NULL ;

    int64_t n = A->n, top, i, j, lnz = 0, jnew, k;
    size_t size;

    int64_t* post = NULL;
    int64_t* c = NULL;

    
    // h is the history vector utilized for the sparse REF
    // triangular solve algorithm. h serves as a global
    // vector which is repeatedly passed into the triangular
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
    // initialize workspace history array
    for (i = 0; i < n; i++)
    {
        h[i] = -1;
    }

    // Obtain the elimination tree of A
    S->parent = IP_Chol_etree(A);         // Obtain the elim tree
    post = IP_Chol_post(S->parent, n);    // Postorder the tree
    
    // Get the column counts of A
    c = IP_Chol_counts(A, S->parent, post, 0);
    
    S->cp = (int64_t*) SLIP_malloc( (n+1)*sizeof(int64_t*));
    S->lnz = IP_cumsum_chol(S->cp, c, n);    // Get column pointers for L
   
    //--------------------------------------------------------------------------
    // allocate and initialize the workspace x
    //--------------------------------------------------------------------------

    // IP_Chol utilizes arbitrary sized integers which can grow beyond the
    // default 64 bits allocated by GMP. If the integers frequently grow, GMP
    // can get bogged down by performing intermediate reallocations. Instead,
    // we utilize a larger estimate on the workspace x vector so that computing
    // the values in L and U do not require too many extra intermediate calls to
    // realloc.
    //
    // Note that the estimate presented here is not an upper bound nor a lower
    // bound.  It is still possible that more bits will be required which is
    // correctly handled internally.
    int64_t estimate = 64 * SLIP_MAX (2, ceil (log2 ((double) n))) ;

    // Create x, a global dense mpz_t matrix of dimension n*1. Unlike rhos, the
    // second boolean parameter is set to false to avoid initializing
    // each mpz entry of x with default size.  It is initialized below.
    SLIP_CHECK (SLIP_matrix_allocate(&x, SLIP_DENSE, SLIP_MPZ, n, 1, n,
        false, /* do not initialize the entries of x: */ false, option));
    
    // Create rhos, a global dense mpz_t matrix of dimension n*1. 
    SLIP_CHECK (SLIP_matrix_allocate(&rhos, SLIP_DENSE, SLIP_MPZ, n, 1, n,
        false, true, option));
    
    if (!x || !rhos)
    {
        FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    
    // initialize the entries of x
    for (i = 0; i < n; i++)
    {
        // Allocate memory for entries of x
        SLIP_CHECK(SLIP_mpz_init2(x->x.mpz[i], estimate));
    }

    //--------------------------------------------------------------------------
    // Declare memory for L 
    //--------------------------------------------------------------------------
    
    // If we are performing a left-looking factorization, we pre-allocate L
    // by performing a symbolic version of the factorization and obtaining the 
    // exact nonzero pattern of L.
    // Conversely, if we are performing an up-looking factorization, we allocate
    // L without initializing eah entry.
    // In both cases, the inidividual (x) values of L are not allocated. Instead,
    // a more efficient method to allocate these values is done inside the 
    // factorization to reduce memory usage.
    
    if (left)
    {
        OK(IP_Pre_Left_Factor(A, &L, xi, S->parent, S, c));
    }
    else
    {
        SLIP_CHECK (SLIP_matrix_allocate(&L, SLIP_CSC, SLIP_MPZ, n, n, S->lnz,
            false, false, option));
    }
    
    // Set the column pointers of L
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
    

    //--------------------------------------------------------------------------
    // Perform the factorization
    //--------------------------------------------------------------------------
    if (left)
    {
        //--------------------------------------------------------------------------
            // Iterations 0:n-1 (1:n in standard)
        //--------------------------------------------------------------------------
        for (k = 0; k < n; k++)
        {
            // LDx = A(:,k)
            OK (IP_Left_Chol_triangular_solve(&top, L, A, k, xi, rhos, h, x, S->parent, c));

            // Set the pivot element
            if (mpz_sgn(x->x.mpz[k]) != 0)
            {
                //pivot = k;
                OK(SLIP_mpz_set(rhos->x.mpz[k], x->x.mpz[k]));
            }
            else
            {
                FREE_WORKSPACE;
                return SLIP_SINGULAR;
            }
            //----------------------------------------------------------------------
            // Iterate accross the nonzeros in x
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
        }
    }
    else
    {
        //--------------------------------------------------------------------------
        // Iterations 0:n-1 (1:n in standard)
        //--------------------------------------------------------------------------
        for (k = 0; k < n; k++)
        {
            // LDx = A(:,k)
            SLIP_CHECK ( IP_Up_Chol_triangular_solve(&top, L, A, k, xi, S->parent, c, rhos, h, x));
        
            // If x[k] is nonzero that is the pivot. if x[k] == 0 then matrix is singular.
            if (mpz_sgn(x->x.mpz[k]) != 0)
            {
                //pivot = k;
                OK(SLIP_mpz_set(rhos->x.mpz[k], x->x.mpz[k]));
            }
            else
            {
                FREE_WORKSPACE;
                return SLIP_SINGULAR;
            }
        
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
