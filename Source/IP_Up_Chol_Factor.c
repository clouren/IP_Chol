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
    SLIP_sparse* A,           // matrix to be factored
    SLIP_sparse* L,           // lower triangular matrix
    Sym_chol* S,           // Permutation, elimination tree etc
    mpz_t* rhos,           // sequence of pivots
    SLIP_options* option// command options
)
{
    SLIP_info ok;
    // Input check
    if (!A || !L || !S || !rhos || !option)
    {
        return SLIP_INCORRECT_INPUT;
    }
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
    mpz_t sigma; SLIP_MPZ_SET_NULL(sigma);
    mpfr_t temp; SLIP_MPFR_SET_NULL(temp);
    mpz_t* x = NULL;
    
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
    // Get most dense column and max of A 
    //--------------------------------------------------------------------------
    // Initialize sigma
    OK(SLIP_mpz_init(sigma)); 
    OK(SLIP_mpz_set(sigma,A->x[0]));
    
    // Get sigma = max(A)
    for (i = 1; i < A->nz; i++)
    {
        if(mpz_cmpabs(sigma,A->x[i]) < 0)
        {
            OK(SLIP_mpz_set(sigma,A->x[i]));
        }
    }
    // sigma = |sigma|
    mpz_abs(sigma,sigma);

    int gamma = A->p[1];
    // get gamma as most dense column
    for (i = 1; i<n; i++)
    {
        if( gamma < A->p[i+1] - A->p[i])
            gamma = A->p[i+1]-A->p[i];
    }
    OK(SLIP_mpfr_init2(temp, 256)); 
    // temp = sigma
    OK(SLIP_mpfr_set_z(temp, sigma, MPFR_RNDN));

    //--------------------------------------------------------------------------
    // Bound = gamma*log2(sigma sqrt(gamma))
    //--------------------------------------------------------------------------
    // temp = sigma*sqrt(gamma)
    OK(SLIP_mpfr_mul_d(temp, temp, (double) sqrt(gamma), MPFR_RNDN));
    // temp = log2(temp)
    OK(SLIP_mpfr_log2(temp, temp, MPFR_RNDN));
    // inner2 = temp
    double inner2 = mpfr_get_d(temp, MPFR_RNDN);
    // Free cache from log2
    mpfr_free_cache();
    // bound = gamma * inner2+1
    int bound = (int) ceil(gamma*(inner2+1));
    // Ensure bound is at least 64 bit
    if (bound < 64) bound = 64;    
    // Free memory
    
    

    //--------------------------------------------------------------------------
    // Declare memory for x, L 
    //--------------------------------------------------------------------------
    // Initialize x
    x = IP_create_mpz_array2(n,bound);
    if (!x)
    {
        FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    
    // Allocate L  
    OK(IP_sparse_alloc2(L, n, n, S->lnz));
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
    
    //--------------------------------------------------------------------------
    // Iteration 0, Set L(0,0)
    //--------------------------------------------------------------------------
    // x = A(:,col)
    OK(IP_Up_get_column(A,0,x)); 
    if (mpz_sgn(x[0]) != 0)
    {
        pivot = 0;
        SLIP_mpz_set(rhos[0], x[0]);
    }
    else
    {
        pivot = SLIP_SINGULAR;
        FREE_WORKSPACE;
        return pivot;
    }
    // top: nnz in column col
    L->i[0] = 0;
    size = mpz_sizeinbase(x[0],2);
    // GMP manual: Allocated size should be size+2
    OK(SLIP_mpz_init2(L->x[0],size+2));
    // Set L[x]
    OK(SLIP_mpz_set(L->x[0],x[0]));
    c[0]++;

    //--------------------------------------------------------------------------
    // Iterations 1:n-1 (2:n in standard)
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
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
            size = mpz_sizeinbase(x[jnew],2);
            // GMP manual: Allocated size should be size+2
            OK(SLIP_mpz_init2(L->x[p], size+2));
            // Place the x value of the L->nz nonzero
            OK(SLIP_mpz_set(L->x[p],x[jnew]));
        }
        p = c[k]++;
        L->i[p] = k;
        size = mpz_sizeinbase(x[k], 2);
        OK(SLIP_mpz_init2(L->x[p], size+2));
        OK(SLIP_mpz_set(L->x[p], x[k]));
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