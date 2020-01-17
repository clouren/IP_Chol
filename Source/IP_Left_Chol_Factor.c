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
    SLIP_sparse* A,           // matrix to be factored
    SLIP_sparse* L,           // lower triangular matrix
    Sym_chol * S,           // stores guess on nnz and column permutation
    mpz_t* rhos,           // sequence of pivots
    int* pinv,             // inverse row permutation
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

    int32_t* h = NULL;
    int32_t* xi = NULL;
    int32_t* c = NULL;
    int32_t* post = NULL;
    mpz_t* x = NULL;
    mpfr_t temp; SLIP_MPFR_SET_NULL(temp);
    mpz_t sigma; SLIP_MPZ_SET_NULL(sigma);
    
    
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
    // Get most dense column and max of A 
    //--------------------------------------------------------------------------
    // Initialize sigma
    mpz_init(sigma); 
    OK(SLIP_mpz_set(sigma,A->x[0]));
    if (ok != SLIP_OK) return ok;
    
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
    SLIP_mpfr_init2(temp, 256); 
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
    for (i = 0; i < n; i++) pinv[i] = i;
    
    // Allocate L  
    OK(IP_Pre_Left_Factor(A, L, xi, S->parent, S, c));
    
    //--------------------------------------------------------------------------
    // Iteration 0, must select pivot
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        c[k] = L->p[k];
    }
    OK(IP_get_column(x, A, 0)); 
    if (mpz_sgn(x[0]) != 0)
    {
        pivot = 0;
        OK(SLIP_mpz_set(rhos[0], x[0]));
    }
    else
        pivot = SLIP_SINGULAR;
    // top: nnz in column col
    top = n-A->p[1]; j = 0;

    // Populate nonzero pattern
    for (i = A->p[0]; i < A->p[1]; i++)
    {
        xi[top+j] = A->i[i];
        j+=1;
    }
    
    // Some error
    if (pivot < SLIP_OK)
    {
        FREE_WORKSPACE;
        return pivot;
    }
    qsort(&xi[top], n-top, sizeof(int32_t), compare3); 
    // Populate L and U
    for (j = top; j < n; j++)
    {
        jnew = xi[j];
        if (jnew >= 0)
        {
            // ith value of x[j]
            size = mpz_sizeinbase(x[jnew],2);
            // GMP manual: Allocated size should be size+2
            OK(SLIP_mpz_init2(L->x[lnz],size+2));
            // Set L[x]
            OK(SLIP_mpz_set(L->x[lnz],x[jnew]));
            lnz += 1;
        }
    }

  
    //--------------------------------------------------------------------------
    // Iterations 1:n-1 (2:n in standard)
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
    {
        // LDx = A(:,k)
        top = IP_Left_Chol_triangular_solve(L, A, k, xi, rhos, pinv, h, x, S->parent, c);
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
                size = mpz_sizeinbase(x[jnew],2);
                // GMP manual: Allocated size should be size+2
                OK(SLIP_mpz_init2(L->x[lnz], size+2));
                // Place the x value of the L->nz nonzero
                OK(SLIP_mpz_set(L->x[lnz],x[jnew]));
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