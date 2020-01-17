#ifndef IP_Chol
#define IP_Chol

#include "../SLIP_LU/v0.6/Headers/SLIP_LU_config.h"

/* This header contains functions to perform the two integer-preserving Cholesky factorizations */

#define SLIP_MARKED(Ap,j) (Ap [j] < 0)
#define SLIP_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SLIP_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SLIP_FLIP(i) (-(i)-2)
#define SLIP_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define SLIP_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
# define SLIP_FREE(x)       \
{                           \
    if ( x != NULL)         \
    {                       \
        SLIP_free (x) ;     \
        x = NULL ;          \
    }                       \
}                           \


typedef struct Sym_chol
{
    int* pinv;      // Row permutation
    int* parent;    // Elimination tree for Cholesky
    int* cp;        // Column pointers for Cholesky
    int lnz;        // Number of nonzeros in Cholesky L
} Sym_chol;
    

    
/* This function frees the Sym_chol data structure */
void Sym_chol_free
(
    Sym_chol* S
)
{
    //SLIP_FREE(S->pinv);
    SLIP_FREE(S->parent);
    SLIP_FREE(S->cp);
    SLIP_FREE(S);
}

   
/* Purpose: This function obtains A(1:k-1,k) and stores it in a
 * dense vector x
 */
void Up_get_column  
(
    SLIP_mat* A,    // input matrix
    int k,          // column to extract
    mpz_t* x        // A(1:k-1,k)
)
{
    // Iterating accross the nonzeros in column k
    for (int i = A->p[k]; i < A->p[k+1]; i++)
    {
        if (A->i[i] <= k)
        {
            // Value of the ith nonzero    
            slip_mpz_set(x[A->i[i]], A->x[i]);
        }
    }
}

/* Permute the matrix A and return A2 = PAP */
SLIP_mat* Chol_permute_A
(
    SLIP_mat* A,        // Initial input matrix
    int* pinv,          // Row permutation
    SLIP_col* S         // Column permutation
)
{
    int nz = 0, j, n = A->n;
    SLIP_mat* A2 = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
    SLIP_mat_alloc(n, n, A->nz, A2);
    for (int k = 0 ; k < n ; k++)
    {
        A2->p [k] = nz ;                       // column k of A2 is column q[k] of A 
        j = S->q [k];
        for (int t = A->p [j] ; t < A->p [j+1] ; t++)
        {
            slip_mpz_set(A2->x[nz], A->x[t]);  // row i of A is row pinv[i] of C 
            A2->i [nz++] = pinv [A->i [t]];
        }
    }
    A2->p [n] = nz ;                       // finalize the last column of C 
    A2->nz = A2->nzmax;
    return A2;
}


/* Compute the elimination tree of A */

int* Chol_etree 
(
    SLIP_mat* A // Input matrix (must be SPD)
)
{
    int i, k, p, m, n, inext, *w, *parent, *ancestor, *prev ;
    if (!A) return (NULL) ;        /* check inputs */
    m = A->m ; n = A->n ;
    parent = (int*) SLIP_malloc(n, sizeof(int));
    w = (int*) SLIP_malloc(n+m, sizeof(int));
    ancestor = w ; prev = w + n ;
    for (k = 0 ; k < n ; k++)
    {
        parent [k] = -1 ;                   /* node k has no parent yet */
        ancestor [k] = -1 ;                 /* nor does k have an ancestor */
        for (p = A->p [k] ; p < A->p [k+1] ; p++)
        {
            i = A->i [p] ;
            for ( ; i != -1 && i < k ; i = inext)   /* traverse from i to k */
            {
                inext = ancestor [i] ;              /* inext = ancestor of i */
                ancestor [i] = k ;                  /* path compression */
                if (inext == -1) parent [i] = k ;   /* no anc., parent is k */
            }
        }
    }
    SLIP_free(w);
    return parent;
}

/* This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. This is more efficient than the SLIP_reach function 
   It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
   
int Chol_ereach 
(
    SLIP_mat *A,    // Matrix to be analyzed
    int k,          // Node to start at
    int* parent,    // ELimination Tree
    int* s,         // Contains the nonzero pattern in s[top..n-1]
    int* w          // Workspace array
)
{
    int i, p, n, len, top ;
    if (!A || !parent || !s || !w) return (-1) ;   /* check inputs */
    top = n = A->n ; 
    SLIP_MARK (w, k) ;                /* mark node k as visited */
    for (p = A->p [k] ; p < A->p [k+1] ; p++)
    {
        i = A->i [p] ;                /* A(i,k) is nonzero */
        if (i > k) continue ;        /* only use upper triangular part of A */
        for (len = 0 ; !SLIP_MARKED (w,i) ; i = parent [i]) /* traverse up etree*/
        {
            s [len++] = i ;         /* L(k,i) is nonzero */
            SLIP_MARK (w, i) ;        /* mark i as visited */
        }
        while (len > 0) s [--top] = s [--len] ; /* push path onto stack */
    }
    for (p = top ; p < n ; p++) SLIP_MARK (w, s [p]) ;    /* unmark all nodes */
    SLIP_MARK (w, k) ;                /* unmark node k */
    return (top) ;                  /* s [top..n-1] contains pattern of L(k,:)*/
}

/* Depth-first search and postorder of a tree rooted at node j */

int Chol_tdfs 
(
    int j,      // Root node
    int k,      
    int* head,  // Head of list
    int* next,  // Next node in the list
    int* post,  // Post ordered tree
    int* stack  // Stack of nodes
)
{ 
    int i, p, top = 0 ;
    if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
    stack [0] = j ;                 /* place j on the stack */
    while (top >= 0)                /* while (stack is not empty) */
    {
        p = stack [top] ;           /* p = top of stack */
        i = head [p] ;              /* i = youngest child of p */
        if (i == -1)
        {
            top-- ;                 /* p has no unordered children left */
            post [k++] = p ;        /* node p is the kth postordered node */
        }
        else
        {
            head [p] = next [i] ;   /* remove i from children of p */
            stack [++top] = i ;     /* start dfs on child node i */
        }
    }
    return (k) ;
}

/* post order a forest */
int *Chol_post 
(
    int* parent,    // Parent[j] is parent of node j in forest
    int n           // Number of nodes in the forest
)
{
    int j, k = 0, *post, *w, *head, *next, *stack ;
    if (!parent) return (NULL) ;                               /* check inputs */
    post = (int*) SLIP_malloc(n, sizeof(int));                 /* allocate result */
    w = (int*) SLIP_malloc (3*n, sizeof (int)) ;               /* get workspace */
    if (!w || !post) return (NULL) ;
    head = w ; next = w + n ; stack = w + 2*n ;
    for (j = 0 ; j < n ; j++) head [j] = -1 ;           /* empty linked lists */
    for (j = n-1 ; j >= 0 ; j--)            /* traverse nodes in reverse order*/
    {
        if (parent [j] == -1) continue ;    /* j is a root */
        next [j] = head [parent [j]] ;      /* add j to list of its parent */
        head [parent [j]] = j ;
    }
    for (j = 0 ; j < n ; j++)
    {
        if (parent [j] != -1) continue ;    /* skip j if it is not a root */
        k = Chol_tdfs (j, k, head, next, post, stack) ;
    }
    SLIP_free(w);
    return (post) ;  /* success; free w, return post */
}


/* consider A(i,j), node j in ith row subtree and return lca(jprev,j) 
   Used to determine Column counts of cholesky factor*/
int Chol_leaf 
(
    int i, 
    int j, 
    int* first, 
    int* maxfirst, 
    int* prevleaf,
    int* ancestor, 
    int* jleaf
)
{
    int q, s, sparent, jprev ;
    if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1) ;
    *jleaf = 0 ;
    if (i <= j || first [j] <= maxfirst [i]) return (-1) ;  /* j not a leaf */
    maxfirst [i] = first [j] ;      /* update max first[j] seen so far */
    jprev = prevleaf [i] ;          /* jprev = previous leaf of ith subtree */
    prevleaf [i] = j ;
    *jleaf = (jprev == -1) ? 1: 2 ; /* j is first or subsequent leaf */
    if (*jleaf == 1) return (i) ;   /* if 1st leaf, q = root of ith subtree */
    for (q = jprev ; q != ancestor [q] ; q = ancestor [q]) ;
    for (s = jprev ; s != q ; s = sparent)
    {
        sparent = ancestor [s] ;    /* path compression */
        ancestor [s] = q ;
    }
    return (q) ;                    /* q = least common ancester (jprev,j) */
}

#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
static void init_ata 
(
    SLIP_mat *AT, 
    int* post, 
    int *w, 
    int **head, 
    int **next
)
{
    int i, k, p, m = AT->n, n = AT->m;
    *head = w+4*n, *next = w+5*n+1 ;
    for (k = 0 ; k < n ; k++) w [post [k]] = k ;    /* invert post */
    for (i = 0 ; i < m ; i++)
    {
        for (k = n, p = AT->p[i] ; p < AT->p[i+1] ; p++) k = CS_MIN (k, w [AT->i[p]]);
        (*next) [i] = (*head) [k] ;     /* place row i in linked list k */
        (*head) [k] = i ;
    }
}

int *Chol_counts 
(
    SLIP_mat *A, 
    int *parent, 
    int *post, 
    int ata // Parameter if we are doing A or A^T A. Setit as 0
)
{
    int i, j, k, n, m, J, s, p, q, jleaf, *maxfirst, *prevleaf,
        *ancestor, *head = NULL, *next = NULL, *colcount, *w, *first, *delta ;
    if (!A || !parent || !post) return (NULL) ;    /* check inputs */
    m = A->m ; n = A->n ;
    s = 4*n + (ata ? (n+m+1) : 0) ;
    delta = colcount = (int*) SLIP_malloc (n, sizeof (int)) ;    /* allocate result */
    w = (int*) SLIP_malloc (s, sizeof (int)) ;                   /* get workspace */
    ancestor = w ; maxfirst = w+n ; prevleaf = w+2*n ; first = w+3*n ;
    for (k = 0 ; k < s ; k++) w [k] = -1 ;      /* clear workspace w [0..s-1] */
    for (k = 0 ; k < n ; k++)                   /* find first [j] */
    {
        j = post [k] ;
        delta [j] = (first [j] == -1) ? 1 : 0 ;  /* delta[j]=1 if j is a leaf */
        for ( ; j != -1 && first [j] == -1 ; j = parent [j]) first [j] = k ;
    }
    for (i = 0 ; i < n ; i++) ancestor [i] = i ; /* each node in its own set */
    for (k = 0 ; k < n ; k++)
    {
        j = post [k] ;          /* j is the kth node in postordered etree */
        if (parent [j] != -1) delta [parent [j]]-- ;    /* j is not a root */
        for (J = HEAD (k,j) ; J != -1 ; J = NEXT (J))   /* J=j for LL'=A case */
        {
            for (p = A->p [J] ; p < A->p [J+1] ; p++)
            {
                i = A->i [p] ;
                q = Chol_leaf (i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
                if (jleaf >= 1) delta [j]++ ;   /* A(i,j) is in skeleton */
                if (jleaf == 2) delta [q]-- ;   /* account for overlap in q */
            }
        }
        if (parent [j] != -1) ancestor [j] = parent [j] ;
    }
    for (j = 0 ; j < n ; j++)           /* sum up delta's of each child */
    {
        if (parent [j] != -1) colcount [parent [j]] += colcount [j] ;
    }
    SLIP_free(w);
    return colcount;    
} 

/* Purpose: This function performs the symmetric sparse REF triangular solve. for uplooking
 * Cholesky factorization. i.e., 
 * (LD) x = A(1:k-1,k). 
 */
int Up_chol_triangular_solve // performs the sparse REF triangular solve
(
    SLIP_mat* L,              // partial L matrix
    SLIP_mat* A,              // input matrix
    int k,                    // iteration of algorithm
    int* xi,                  // nonzero pattern vector
    int* parent,              // Elimination tree
    int* c,                   // Column pointers
    mpz_t* rhos,              // sequence of pivots
    int* h,                   // history vector
    mpz_t* x                  // solution of system ==> kth column of L and U
)
{
    if (!L || !A || !xi || !parent || !c || !rhos || !h || !x)
        return SLIP_INCORRECT_INPUT;
    int j, jnew, i, inew, p, m, top, n = A->n, col, ok;
    
    //--------------------------------------------------------------------------
    // Initialize REF TS by getting nonzero patern of x && obtaining A(:,k)
    //--------------------------------------------------------------------------
    top = Chol_ereach(A, k, parent, xi, c);  // Obtain nonzero pattern in xi[top..n]
    std::sort(xi + top, xi + n);             // Sort the nonzero pattern. No row perm is applied so can use simple sort
    SLIP_reset_mpz_array_2(x,n,top,xi);      // Reset x[i] = 0 for all i in nonzero pattern
    mpz_set_ui(x[k], 0);                     // Reset x[k]
    SLIP_reset_int_array_2(h,n,top,xi);      // Reset h[i] = -1 for all i in nonzero pattern
    Up_get_column(A,k,x);                    // Set x = A(:,k)
    
    //--------------------------------------------------------------------------
    // Iterate accross nonzeros in x
    //--------------------------------------------------------------------------
    ok = SLIP_OK;                       // TODO: Implement better way to do this
    for (p = top; p < n; p++)
    {   
        if (ok != SLIP_OK) break;
        /* Finalize x[j] */
        j = xi[p];                          // First nonzero term
        //if (j == k) continue;             // Do not operate on x[k] in TS
        if (mpz_sgn(x[j]) == 0) continue;   // If x[j] == 0 no work must be done
                
        // History update x[j]
        if (h[j] < j-1)
        {
            // History update x[j]
            ok = slip_mpz_mul(x[j], x[j], rhos[j-1]);
            if (h[j] > -1)
            {
                ok = slip_mpz_divexact(x[j], x[j], rhos[ h[j]]);
            }
        }
        
        //------------------------------------------------------------------
        // IPGE updates
        //------------------------------------------------------------------
        // ----------- Iterate accross nonzeros in Lij ---------------------
        for (m = L->p[j] +1; m < c[j]; m++)
        {
            i = L->i[m];            // i value of Lij
            if (i > j && i < k)     // Update all dependent x[i] excluding x[k]
            {
                    /*************** If lij==0 then no update******************/
                if (mpz_sgn(L->x[m]) == 0) continue;

                //----------------------------------------------------------
                /************* lij is nonzero, x[i] is zero****************/
                // x[i] = 0 then only perform IPGE update subtraction/division
                //----------------------------------------------------------
                if (mpz_sgn(x[i]) == 0)
                {
                    // No previous pivot
                    if (j < 1)
                    {
                        ok = slip_mpz_submul(x[i],L->x[m],x[j]);// x[i] = 0 - lij*x[j]
                        if (ok != SLIP_OK) break;
                        h[i] = j;                  // Entry is up to date
                    }
                        
                    // Previous pivot exists
                    else
                    {
                        ok = slip_mpz_submul(x[i],L->x[m],x[j]);// x[i] = 0 - lij*x[j]
                        if (ok != SLIP_OK) break;
                        ok = slip_mpz_divexact(x[i],x[i],rhos[j-1]);// x[i] = x[i] / rho[j-1]
                        if (ok != SLIP_OK) break;
                        h[i] = j;                  // Entry is up to date
                    }
                }

                //----------------------------------------------------------
                /************ Both lij and x[i] are nonzero****************/
                // x[i] != 0 --> History & IPGE update on x[i]
                //----------------------------------------------------------
                else
                {
                    // No previous pivot in this case
                    if (j < 1)
                    {
                        ok = slip_mpz_mul(x[i],x[i],rhos[0]);      // x[i] = x[i]*rho[0]
                        if (ok != SLIP_OK) break;
                        ok = slip_mpz_submul(x[i], L->x[m], x[j]);// x[i] = x[i] - lij*xj
                        if (ok != SLIP_OK) break;
                        h[i] = j;                  // Entry is now up to date
                    }
                    // There is a previous pivot
                    else
                    {
                        // History update if necessary
                        if (h[i] < j - 1)
                        {
                            ok = slip_mpz_mul(x[i],x[i],rhos[j-1]);// x[i] = x[i] * rho[j-1]
                            if (ok != SLIP_OK) break;
                            if (h[i] > -1)
                            {
                                ok = slip_mpz_divexact(x[i],x[i],rhos[h[i]]);// x[i] = x[i] / rho[h[i]]
                                if (ok != SLIP_OK) break;
                            }
                        }
                        ok = slip_mpz_mul(x[i],x[i],rhos[j]);// x[i] = x[i] * rho[j]
                        if (ok != SLIP_OK) break;
                        ok = slip_mpz_submul(x[i], L->x[m], x[j]);// x[i] = x[i] - lij*xj
                        if (ok != SLIP_OK) break;
                        ok = slip_mpz_divexact(x[i],x[i],rhos[j-1]);// x[i] = x[i] / rho[j-1] 
                        if (ok != SLIP_OK) break;
                        h[i] = j;                  // Entry is up to date
                    }
                }
            }
        }
        // Update x[k]
        if (h[k] < j - 1)
        {
            ok = slip_mpz_mul(x[k],x[k],rhos[j-1]); // x[k] = x[k] * rho[j-1]
            if (h[k] > -1)
            {
                ok = slip_mpz_divexact(x[k],x[k],rhos[h[k]]);// x[k] = x[k] / rho[h[k]]
                if (ok != SLIP_OK) break;
            }
        }
        ok = slip_mpz_mul(x[k],x[k],rhos[j]); // x[k] = x[k] * rho[j]
        if (ok != SLIP_OK) break;
        ok = slip_mpz_submul(x[k], x[j], x[j]);// x[k] = x[k] - xj*xj
        if (ok != SLIP_OK) break;
        if (j != 0)
            ok = slip_mpz_divexact(x[k],x[k],rhos[j-1]); // x[k] = x[k] / rho[j-1] 
        if (ok != SLIP_OK) break;
        h[k] = j;   
    }
    // Finalize x[k]
    if (h[k] < k-1)
    {
        ok = slip_mpz_mul(x[k], x[k], rhos[k-1]);
        if (h[k] > -1)
        {
            ok = slip_mpz_divexact(x[k], x[k], rhos[ h[k]]);
        }
    }
    if (ok != SLIP_OK)
        return ok;
    else
        return top;
}

/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
int Up_Chol_Factor         // performs the Up looking Cholesky factorization
(
    SLIP_mat* A,           // matrix to be factored
    SLIP_mat* L,           // lower triangular matrix
    Sym_chol* S,           // Permutation, elimination tree etc
    mpz_t* rhos,           // sequence of pivots
    SLIP_LU_Options* option// command options
)
{
    // Input check
    if (!A || !L || !S || !rhos || !option)
        return SLIP_INCORRECT_INPUT;
    //--------------------------------------------------------------------------
    // Declare and initialize workspace 
    //--------------------------------------------------------------------------
    // Begin timing factorization
    std::chrono::steady_clock::time_point t_begin
        = std::chrono::steady_clock::now();
    int n = A->n, ok, numRealloc = 0, k = 0, top, i, j, col, loc, lnz = 0, unz = 0, pivot,
        check, jnew, *xi, *h, *col_loc, *pivs;
    long size;

    // History vector
    h = (int*) SLIP_malloc(n, sizeof(int));
    // Nonzero pattern
    xi = (int*) SLIP_malloc(2*n, sizeof(int));

    if (!h || !xi) return SLIP_OUT_OF_MEMORY;
    SLIP_reset_int_array(h,n);

    // Obtain the elimination tree of A
    S->parent = Chol_etree(A);              // Obtain the elim tree
    int* post = Chol_post(S->parent, n);    // Postorder the tree
    
    // Get the column counts of A
    int* c = Chol_counts(A, S->parent, post, 0);
    
    S->cp = (int*) SLIP_malloc(n+1, sizeof(int));
    S->lnz = cs_cumsum(S->cp, c, n);    // Get column pointers for L
   
    SLIP_FREE(post); 

    //--------------------------------------------------------------------------
    // Get most dense column and max of A 
    //--------------------------------------------------------------------------
    // Initialize sigma
    mpz_t sigma; mpz_init(sigma); 
    ok = slip_mpz_set(sigma,A->x[0]);
    if (ok != SLIP_OK) return ok;
    
    // Get sigma = max(A)
    for (i = 1; i < A->nz; i++)
    {
        if(mpz_cmpabs(sigma,A->x[i]) < 0)
        {
            ok = slip_mpz_set(sigma,A->x[i]);
            if (ok != SLIP_OK) return ok;
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
    mpfr_t temp; mpfr_init2(temp, 256); 
    // temp = sigma
    ok = slip_mpfr_set_z(temp, sigma, MPFR_RNDN);
    if (ok != SLIP_OK) return ok;

    //--------------------------------------------------------------------------
    // Bound = gamma*log2(sigma sqrt(gamma))
    //--------------------------------------------------------------------------
    // temp = sigma*sqrt(gamma)
    ok = slip_mpfr_mul_d(temp, temp, (double) sqrt(gamma), MPFR_RNDN);
    if (ok != SLIP_OK) return ok;
    // temp = log2(temp)
    ok = slip_mpfr_log2(temp, temp, MPFR_RNDN);
    if (ok != SLIP_OK) return ok;
    // inner2 = temp
    double inner2 = mpfr_get_d(temp, MPFR_RNDN);
    // Free cache from log2
    mpfr_free_cache();
    // bound = gamma * inner2+1
    int bound = std::ceil(gamma*(inner2+1));
    // Ensure bound is at least 64 bit
    if (bound < 64) bound = 64;    
    // Free memory
    mpfr_clear(temp); mpz_clear(sigma);
    

    //--------------------------------------------------------------------------
    // Declare memory for x, L 
    //--------------------------------------------------------------------------
    // Initialize x
    mpz_t* x = SLIP_initialize_mpz_array2(n,bound);
    if (!x) return SLIP_OUT_OF_MEMORY;
    
    // Allocate L  
    ok = SLIP_mat_alloc2(n, n, S->lnz, L);
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
    if (ok != SLIP_OK) return ok;
    
    //--------------------------------------------------------------------------
    // Iteration 0, Set L(0,0)
    //--------------------------------------------------------------------------
    // x = A(:,col)
    Up_get_column(A,0,x); 
    if (mpz_sgn(x[0]) != 0)
    {
        pivot = 0;
        slip_mpz_set(rhos[0], x[0]);
    }
    else
    {
        pivot = SLIP_SINGULAR;
    }
    // top: nnz in column col
    L->i[0] = 0;
    size = mpz_sizeinbase(x[0],2);
    // GMP manual: Allocated size should be size+2
    ok = slip_mpz_init2(L->x[0],size+2);
    // Set L[x]
    ok = slip_mpz_set(L->x[0],x[0]);
    c[0]++;

    //--------------------------------------------------------------------------
    // Iterations 1:n-1 (2:n in standard)
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
    {
        if (ok != SLIP_OK) break;
                
        // LDx = A(:,k)
        top = Up_chol_triangular_solve(L, A, k, xi, S->parent, c, rhos, h, x);
        if (top < 0)
        {
            ok = top;
            break;
        }
        if (mpz_sgn(x[k]) != 0)
        {
            pivot = k;
            mpz_set(rhos[k], x[k]);
        }
        else
        {
            pivot = SLIP_SINGULAR;
        }
        // Error
        if (pivot < SLIP_OK)
        {
            ok = pivot;
            break;
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
            ok = slip_mpz_init2(L->x[p], size+2);
            if (ok != SLIP_OK) break;
            // Place the x value of the L->nz nonzero
            ok = slip_mpz_set(L->x[p],x[jnew]);
            if (ok != SLIP_OK) break;
        }
        p = c[k]++;
        L->i[p] = k;
        size = mpz_sizeinbase(x[k], 2);
        ok = slip_mpz_init2(L->x[p], size+2);
        ok = slip_mpz_set(L->x[p], x[k]);
    }
    L->nz = S->lnz;
    // Finalize L->p
    L->p[n] = L->nz;        

    //--------------------------------------------------------------------------
    // find the elasped time
    //--------------------------------------------------------------------------
    std::chrono::steady_clock::time_point t_end 
        = std::chrono::steady_clock::now();
    option->t_f = std::chrono::duration_cast<std::chrono::duration<float>>
        (t_end - t_begin);

    //--------------------------------------------------------------------------
    // Set the determinant of A (may be scaled)
    //--------------------------------------------------------------------------
    if (ok == SLIP_OK)
        ok = slip_mpz_set(option->determinant, rhos[n-1]);
    // Print statistics if desired
    if (option->print == 1)
    {
        std::cout<<"\n\n****Factorization Statistics****";
        std::cout<<"\n\nNum Nonzeros in A: "<<A->nz;
        std::cout<<"\nDimension of A, L: "  <<n;
        std::cout<<"\nEstimated L nonzeros: "<<S->lnz;
        std::cout<<"\nActual L nonzeros: "<<L->nz;
        std::cout<<"\nNum Costly Reallocs: "<<numRealloc;
        std::cout<<"\n\n";
    }

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    SLIP_delete_mpz_array(x,n);
    SLIP_FREE(c);    
    SLIP_free(xi); SLIP_free(h); 
    // Collapse L
    if (ok == SLIP_OK)
        SLIP_mat_collapse(L); 
    return ok;
}

/* This solves the system L'x = b for Cholesky factorization */
int Chol_ltsolve 
(
    SLIP_mat *L,    // The lower triangular matrix
    mpz_t **x,      // Solution vector
    int numRHS      // Number of RHS vectors
)
{
    int p, j, n, k;
    n = L->n;
    for (int k = 0; k < numRHS; k++)
    {
        for (j = n-1; j >= 0; j--)
        {
            if (mpz_sgn(x[j][k]) == 0) continue;
            for (p = L->p[j]+1; p < L->p[j+1]; p++)
            {
                slip_mpz_submul(x[j][k], L->x[p], x[L->i[p]][k]);
            }
            slip_mpz_divexact( x[j][k], x[j][k], L->x[L->p[j]]);
        }
    }
    return SLIP_OK;
}

/* Purpose: This function solves the linear system LD^(-1)L' x = b.*/
int Up_Solve               //solves the linear system LD^(-1)L' x = b
(
    mpq_t** x,              // rational solution to the system
    mpz_t** b,              // right hand side vector
    mpz_t* rhos,            // sequence of pivots
    SLIP_mat* L,            // lower triangular matrix
    int* pinv,              // row permutation
    SLIP_LU_Options* option,// command options
    int numRHS              // number of RHS vectors
)
{
    if (!x || !b || !rhos || !pinv) return SLIP_INCORRECT_INPUT;
    std::chrono::steady_clock::time_point t_s_begin
        = std::chrono::steady_clock::now();
    int i, k, n = L->n, ok = SLIP_OK;
    // Permuted b
    mpz_t** b2 = SLIP_initialize_mpz_mat(n, numRHS);
    if (!b2) return SLIP_OUT_OF_MEMORY;
    // Set workspace b
    for (k = 0; k < numRHS; k++)
    {
        if (ok != SLIP_OK) break;
        for (i = 0; i < n; i++)   
        {
            ok = slip_mpz_set(b2[pinv[i]][k], b[i][k]);
            if (ok != SLIP_OK) break;
        }
    }
    // L*b2 = b2
    if (ok == SLIP_OK);
        ok = SLIP_Forward_Sub(L, b2, rhos, numRHS);
    // b2 = b2 * det 
    if (ok == SLIP_OK);
        ok = SLIP_array_mul(b2, rhos[n-1], n, numRHS);
    // U b2 = b2
    if (ok == SLIP_OK);
        ok = Chol_ltsolve(L, b2, numRHS);
    // x = b2/det 
    if (ok == SLIP_OK);
        ok = SLIP_array_div(x, b2, rhos[n-1], n, numRHS);
    std::chrono::steady_clock::time_point t_s_end
        = std::chrono::steady_clock::now();
    option->t_s = std::chrono::duration_cast<std::chrono::duration<float>>
        (t_s_end - t_s_begin);
    // Free memory
    SLIP_delete_mpz_mat(b2, n, numRHS);
    return ok;
}

/* ========================================================================== */
/* This function prints statistics about the Up-Cholesky factorization            */
/* ========================================================================== */
int Chol_Print_Stats     // prints statistics about SLIP LU factorization
(
    SLIP_LU_Options* option, // option struct telling how much info to print
    int check,               // whether the solution is correct or not
    int nnz,                 // number of nonzeros in L+U
    int n,                   // dimension of A
    int numRHS,              // number of RHS vectors
    mpq_t** x                // final solution vector
)
{
    if (!option || !x) return SLIP_INCORRECT_INPUT;
    if (option->status != SLIP_OK) return option->status;
    // Was the solution checked?
    if (option->check == 1)
    {
        if (check == SLIP_CORRECT)
            std::cout<<"\nSolution is correct!";

        // Incorrect solution. This should not happen
        else
        {
            std::cout<<"\n****ERROR**** Solution incorrect!";
            std::cout<<"\n\nHave you modified the source code?";
            std::cout<<"\nReinstall and if this issue persists email me at:"
                <<" clouren@tamu.edu\n\n";
        }
        std::cout<<"\nSolution check time:\t\t"<<option->t_c.count();
    }

    // Total Run time
    option->t_tot = option->t_f + option->t_s+option->t_inp;
    
    // Info about output file
    if (option->print2 == 1)
    {
        std::cout<<"\nSolution output to file named:\t"<< option->out_file;
        if (option->rat == 1)
            std::cout<<"\nSolution output in "
            <<"full precision rational arithmetic";
        else if (option->rat == 2)
            std::cout<<"\nSolution output in double precision";
        else if (option->rat == 3)
            std::cout<<"\nSolution output in fixed precision of size: "
            << option->prec << " bits";
        SLIP_LU_Print_File(option, x, n, numRHS);
    }

    // Print out factorization statistics
    if (option->print3 == 1)
    {
        double nsquare = (double) n*n;
        std::cout<<"\nNumber of L nonzeros:\t\t"       << nnz;
        std::cout<<"\nL Density:\t\t\t"                << (double)nnz/(nsquare);
        std::cout<<"\nInput manipulation time:\t"      << option->t_inp.count();
        std::cout<<"\nSymbolic Analysis time:\t\t"     << option->t_sym.count();
        std::cout<<"\nUp-Chol Factorization time:\t"   << option->t_f.count();
        std::cout<<"\nUp-Chol F/B Substitution time:\t"<< option->t_s.count();
        std::cout<<"\nUp-Chol Total Run Time:\t\t"     << option->t_tot.count()
            <<"\n\n";
    }
    return SLIP_OK;
}

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c 
   From Tim Davis SuiteSparse */
double cs_cumsum_chol 
(
    int *p, 
    int *c, 
    int n
)
{
    int i, nz = 0 ;
    double nz2 = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        nz2 += c [i] ;              /* also in double to avoid int overflow */
        c [i] = p [i] ;             /* also copy p[0..n-1] back into c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz2) ;                  /* return sum (c [0..n-1]) */
}

/* This function sets C = A' */
SLIP_mat* SLIP_transpose_chol 
(
    SLIP_mat *A     // Matrix to be transposed
)
{
    SLIP_mat* C = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
    SLIP_mat_alloc (A->n, A->m, A->nz, C);
    int p, q, j, n, m, *w ;
    m = A->m ; n = A->n ; 
    w = (int*) SLIP_malloc(m, sizeof(int));
    for (p = 0; p < m; p++) w[p] = 0;
    for (p = 0 ; p < A->p [n] ; p++) w [A->i [p]]++ ;       /* row counts */
    cs_cumsum_chol (C->p, w, m) ;                               /* row pointers */
    for (j = 0 ; j < n ; j++)
    {
        for (p = A->p [j] ; p < A->p [j+1] ; p++)
        {
            q = w [A->i [p]]++;
            C->i [q] = j ;                 /* place A(i,j) as entry C(j,i) */
            mpz_set(C->x[q], A->x[p]);
        }
    }
    C->nz = A->nz; 
    C->p[m] = C->nz;
    SLIP_free(w);
    return C;
}

/* Purpose: This function solves the linear system LD^(-1)L' x = b.*/
int SLIP_Chol_Solve             //solves the linear system LD^(-1)L' x = b
(
    mpq_t** x,              // rational solution to the system
    mpz_t** b,              // right hand side vector
    mpz_t* rhos,            // sequence of pivots
    SLIP_mat* L,            // lower triangular matrix
    SLIP_mat* U,            // upper triangular matrix
    int* pinv,              // row permutation
    SLIP_LU_Options* option,// command options
    int numRHS              // number of RHS vectors
)
{
    if (!x || !b || !rhos || !pinv) return SLIP_INCORRECT_INPUT;
    std::chrono::steady_clock::time_point t_s_begin
        = std::chrono::steady_clock::now();
    int i, k, n = L->n, ok = SLIP_OK;
    // Permuted b
    mpz_t** b2 = SLIP_initialize_mpz_mat(n, numRHS);
    if (!b2) return SLIP_OUT_OF_MEMORY;
    // Set workspace b
    for (k = 0; k < numRHS; k++)
    {
        if (ok != SLIP_OK) break;
        for (i = 0; i < n; i++)   
        {
            ok = slip_mpz_set(b2[pinv[i]][k], b[i][k]);
            if (ok != SLIP_OK) break;
        }
    }
    // L*b2 = b2
    if (ok == SLIP_OK);
        ok = SLIP_Forward_Sub(L, b2, rhos, numRHS);
    // b2 = b2 * det 
    if (ok == SLIP_OK);
        ok = SLIP_array_mul(b2, rhos[n-1], n, numRHS);
    // U b2 = b2
    if (ok == SLIP_OK);
        ok = SLIP_Back_Sub(U, b2, numRHS);
    // x = b2/det 
    if (ok == SLIP_OK);
        ok = SLIP_array_div(x, b2, rhos[n-1], n, numRHS);
    std::chrono::steady_clock::time_point t_s_end
        = std::chrono::steady_clock::now();
    option->t_s = std::chrono::duration_cast<std::chrono::duration<float>>
        (t_s_end - t_s_begin);
    // Free memory
    SLIP_delete_mpz_mat(b2, n, numRHS);
    return ok;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//--------------------------Alternate Left looking----------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------


/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
int Pre_Left_Factor         // performs the Up looking Cholesky factorization
(
    SLIP_mat* A,
    SLIP_mat* L,              // partial L matrix
    int* xi,                  // nonzero pattern vector
    int* parent,              // Elimination tree
    Sym_chol * S,           // stores guess on nnz and column permutation
    int* c                   // Column pointers
)
{
    // Input check/
    if (!L || !xi || !parent)
        return SLIP_INCORRECT_INPUT;
    
    int top, k, ok, i, j, jnew, n = A->n;
    //--------------------------------------------------------------------------
    // Declare memory for L 
    //--------------------------------------------------------------------------
       
    // Allocate L  
    ok = SLIP_mat_alloc2(n, n, S->lnz, L);
    for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
        
    L->i[0] = 0;
    c[0]++;

    //--------------------------------------------------------------------------
    // Iterations 1:n-1 (2:n in standard)
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
    {
        if (ok != SLIP_OK) break;
        top = Chol_ereach(A, k, parent, xi, c);  // Obtain nonzero pattern in xi[top..n]
        //std::sort(xi + top, xi + n);             // Sort the nonzero pattern. No row perm is applied so can use simple sort
    
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
    return ok;
}


/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
int SLIP_chol_triangular_solve // performs the sparse REF triangular solve
(
    SLIP_mat* L,              // partial L matrix
    SLIP_mat* A,              // input matrix
    int k,                    // iteration of algorithm
    int* xi,                  // nonzero pattern vector
    mpz_t* rhos,              // sequence of pivots
    int* pinv,                // inverse row permutation
    int* h,                   // history vector
    mpz_t* x,                  // solution of system ==> kth column of L and U
    int* parent,
    int* c
)
{
    if (!L || !A || !xi || !rhos || !pinv || !h || !x)
        return SLIP_INCORRECT_INPUT;
    int j, jnew, i, inew, p, m, top, n, col, ok;

    //--------------------------------------------------------------------------
    // Initialize REF TS by getting nonzero patern of x && obtaining A(:,k)
    //--------------------------------------------------------------------------
    n = A->n;                                // Size of matrix and the dense vectors
   
    /* Chol_ereach gives the nonzeros located in L(k,:) upon completion
     * the vector xi contains the indices of the first k-1 nonzeros in column
     * k of L 
     */
    top = Chol_ereach(A, k, parent, xi, c);
    
    j = top; // Store where the first k-1 nonzeros end
    
    // Populate the rest of the nonzero pattern
    for (i = L->p[k]; i < L->p[k+1]; i++)
    {
        top -= 1;           // One more nonzero in column k
        xi[top] = L->i[i];  // Index of the new nonzero
    }
    
    // Reset the array
    SLIP_reset_mpz_array_2(x, n, top, xi);
    
    // Now we obtain the values of the first k-1 entries of x
    for (i = j; i < n; i++)
    {
        m = xi[i];
        p = c[m]++;
        p+=1;
        mpz_set(x[m], L->x[p]);
    }

    for (i = A->p[k]; i < A->p[k+1]; i++)
    {
        if ( A->i[i] >= k)
        {
            mpz_set(x[A->i[i]], A->x[i]);
        }
    }
    std::sort(xi+top, xi+n); 
    SLIP_reset_int_array_2(h,n,top,xi);      // Reset h[i] = -1 for all i in nonzero pattern
        
    //--------------------------------------------------------------------------
    // Iterate accross nonzeros in x
    //--------------------------------------------------------------------------
    ok = SLIP_OK;
    for ( p = top; p < n; p++)
    {   
        if (ok != SLIP_OK) break;
        /* Finalize x[j] */
        j = xi[p];                        // First nonzero term
        if (mpz_sgn(x[j]) == 0) continue;// If x[j] == 0 no work must be done
        if (j < k)                    // jnew < k implies entries in U
        {
            //------------------------------------------------------------------
            // IPGE updates
            //------------------------------------------------------------------
            // ----------- Iterate accross nonzeros in Lij ---------------------
            for (m = L->p[j]; m < L->p[j+1]; m++)
            {
                i = L->i[m];            // i value of Lij
                if (i > j && i >= k)
                {
                    /*************** If lij==0 then no update******************/
                    if (mpz_sgn(L->x[m]) == 0) continue;

                    //----------------------------------------------------------
                    /************* lij is nonzero, x[i] is zero****************/
                    // x[i] = 0 then only perform IPGE update subtraction/division
                    //----------------------------------------------------------
                    if (mpz_sgn(x[i]) == 0)
                    {
                        // No previous pivot
                        if (j < 1)
                        {
                            ok = slip_mpz_submul(x[i],L->x[m],x[j]);// x[i] = 0 - lij*x[j]
                            if (ok != SLIP_OK) break;
                            h[i] = j;                  // Entry is up to date
                        }
                        
                        // Previous pivot exists
                        else
                        {
                            ok = slip_mpz_submul(x[i],L->x[m],x[j]);// x[i] = 0 - lij*x[j]
                            if (ok != SLIP_OK) break;
                            ok = slip_mpz_divexact(x[i],x[i],rhos[j-1]);// x[i] = x[i] / rho[j-1]
                            if (ok != SLIP_OK) break;
                            h[i] = j;                  // Entry is up to date
                        }
                    }

                    //----------------------------------------------------------
                    /************ Both lij and x[i] are nonzero****************/
                    // x[i] != 0 --> History & IPGE update on x[i]
                    //----------------------------------------------------------
                    else
                    {
                        // No previous pivot in this case
                        if (j < 1)
                        {
                            ok = slip_mpz_mul(x[i],x[i],rhos[0]);      // x[i] = x[i]*rho[0]
                            if (ok != SLIP_OK) break;
                            ok = slip_mpz_submul(x[i], L->x[m], x[j]);// x[i] = x[i] - lij*xj
                            if (ok != SLIP_OK) break;
                            h[i] = j;                  // Entry is now up to date
                        }
                        // There is a previous pivot
                        else
                        {
                            // History update if necessary
                            if (h[i] < j - 1)
                            {
                                ok = slip_mpz_mul(x[i],x[i],rhos[j-1]);// x[i] = x[i] * rho[j-1]
                                if (ok != SLIP_OK) break;
                                if (h[i] > -1)
                                {
                                    ok = slip_mpz_divexact(x[i],x[i],rhos[h[i]]);// x[i] = x[i] / rho[h[i]]
                                    if (ok != SLIP_OK) break;
                                }
                            }
                            ok = slip_mpz_mul(x[i],x[i],rhos[j]);// x[i] = x[i] * rho[j]
                            if (ok != SLIP_OK) break;
                            ok = slip_mpz_submul(x[i], L->x[m], x[j]);// x[i] = x[i] - lij*xj
                            if (ok != SLIP_OK) break;
                            ok = slip_mpz_divexact(x[i],x[i],rhos[j-1]);// x[i] = x[i] / rho[j-1] 
                            if (ok != SLIP_OK) break;
                            h[i] = j;                  // Entry is up to date
                        }
                    }
                }
            }
        }
        else                                              // Entries of L
        {
            //------------------------------------------------------------------
            // History update
            //------------------------------------------------------------------
            if (h[j] < k-1)
            {
                ok = slip_mpz_mul(x[j],x[j],rhos[k-1]);           // x[j] = x[j] * rho[k-1]
                if (ok != SLIP_OK) break;
                if (h[j] > -1)
                {
                    ok = slip_mpz_divexact(x[j],x[j],rhos[h[j]]);// x[j] = x[j] / rho[h[j]]
                    if (ok != SLIP_OK) break;
                }
            }
        }
    }
    // Output the beginning of nonzero pattern
    if (ok != SLIP_OK)
        return ok;
    else
        return top;
}

/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
int SLIP_Chol_Factor         // performs the SLIP LU factorization
(
    SLIP_mat* A,           // matrix to be factored
    SLIP_mat* L,           // lower triangular matrix
    Sym_chol * S,           // stores guess on nnz and column permutation
    mpz_t* rhos,           // sequence of pivots
    int* pinv,             // inverse row permutation
    SLIP_LU_Options* option// command options
)
{
    // Input check
    if (!A || !L || !S || !rhos || !pinv || !option)
        return SLIP_INCORRECT_INPUT;
    //--------------------------------------------------------------------------
    // Declare and initialize workspace 
    //--------------------------------------------------------------------------
    // Begin timing factorization
    std::chrono::steady_clock::time_point t_begin
        = std::chrono::steady_clock::now();
    int n, ok, numRealloc = 0, k = 0, top, i, j, col, loc, lnz = 0, unz = 0, pivot,
        check, jnew, *xi, *h, *col_loc, *pivs, *row_perm;
    long size;

    n = A->n;
    // History vector
    h = (int*) SLIP_malloc(n, sizeof(int));
    // Nonzero pattern
    xi = (int*) SLIP_malloc(2*n, sizeof(int));
    // Row permutation, inverse of pinv
    if (!h || !xi) return SLIP_OUT_OF_MEMORY;
    SLIP_reset_int_array(h,n);
    
    // Obtain the elimination tree of A
    S->parent = Chol_etree(A);              // Obtain the elim tree
    int* post = Chol_post(S->parent, n);    // Postorder the tree
    
    // Get the column counts of A
    int* c = Chol_counts(A, S->parent, post, 0);
    
    S->cp = (int*) SLIP_malloc(n+1, sizeof(int));
    S->lnz = cs_cumsum(S->cp, c, n);    // Get column pointers for L
   
    SLIP_FREE(post); 

    //--------------------------------------------------------------------------
    // Get most dense column and max of A 
    //--------------------------------------------------------------------------
    // Initialize sigma
    mpz_t sigma; mpz_init(sigma); 
    ok = slip_mpz_set(sigma,A->x[0]);
    if (ok != SLIP_OK) return ok;
    
    // Get sigma = max(A)
    for (i = 1; i < A->nz; i++)
    {
        if(mpz_cmpabs(sigma,A->x[i]) < 0)
        {
            ok = slip_mpz_set(sigma,A->x[i]);
            if (ok != SLIP_OK) return ok;
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
    mpfr_t temp; mpfr_init2(temp, 256); 
    // temp = sigma
    ok = slip_mpfr_set_z(temp, sigma, MPFR_RNDN);
    if (ok != SLIP_OK) return ok;

    //--------------------------------------------------------------------------
    // Bound = gamma*log2(sigma sqrt(gamma))
    //--------------------------------------------------------------------------
    // temp = sigma*sqrt(gamma)
    ok = slip_mpfr_mul_d(temp, temp, (double) sqrt(gamma), MPFR_RNDN);
    if (ok != SLIP_OK) return ok;
    // temp = log2(temp)
    ok = slip_mpfr_log2(temp, temp, MPFR_RNDN);
    if (ok != SLIP_OK) return ok;
    // inner2 = temp
    double inner2 = mpfr_get_d(temp, MPFR_RNDN);
    // Free cache from log2
    mpfr_free_cache();
    // bound = gamma * inner2+1
    int bound = std::ceil(gamma*(inner2+1));
    // Ensure bound is at least 64 bit
    if (bound < 64) bound = 64;    
    // Free memory
    mpfr_clear(temp); mpz_clear(sigma);
    

    //--------------------------------------------------------------------------
    // Declare memory for x, L 
    //--------------------------------------------------------------------------
    // Initialize x
    mpz_t* x = SLIP_initialize_mpz_array2(n,bound);
    if (!x) return SLIP_OUT_OF_MEMORY;
    for (i = 0; i < n; i++) pinv[i] = i;
    
    // Allocate L  
    ok = Pre_Left_Factor(A, L, xi, S->parent, S, c);
    if (ok != SLIP_OK) return ok;
    
    //--------------------------------------------------------------------------
    // Iteration 0, must select pivot
    //--------------------------------------------------------------------------
    check = 0;
    for (k = 0; k < n; k++)
    {
        c[k] = L->p[k];
    }
    //std::cout<<"\nc[0] is: "<<c[0];
    // x = A(:,col)
    SLIP_get_column(A,0,x); 
    if (mpz_sgn(x[0]) != 0)
    {
        pivot = 0;
        slip_mpz_set(rhos[0], x[0]);
    }
    else
        pivot = SLIP_SINGULAR;
    // top: nnz in column col
    top = n-A->p[1]; j = 0;

    // Populate nonzero pattern
    for (i = A->p[0]; i < A->p[1]; i++)
    {
        xi[top+j] = A->i[i];
        //if (A->i[i] == 0)
        //{
          //  c[A->i[i]]++;
        //}
        j+=1;
    }
    
    // Some error
    if (pivot < SLIP_OK)
        return pivot;
    std::sort(xi + top, xi + n);
    // Populate L and U
    for (j = top; j < n; j++)
    {
        jnew = xi[j];
        if (jnew >= 0)
        {
            // ith value of x[j]
            size = mpz_sizeinbase(x[jnew],2);
            // GMP manual: Allocated size should be size+2
            ok = slip_mpz_init2(L->x[lnz],size+2);
            if (ok != SLIP_OK) break;
            // Set L[x]
            ok = slip_mpz_set(L->x[lnz],x[jnew]);
            if (ok != SLIP_OK) break;
            lnz += 1;
        }
    }

  
    //--------------------------------------------------------------------------
    // Iterations 1:n-1 (2:n in standard)
    //--------------------------------------------------------------------------
    for (k = 1; k < n; k++)
    {
        if (ok != SLIP_OK) break;
                     
        // LDx = A(:,k)
        top = SLIP_chol_triangular_solve(L, A, k, xi, rhos, pinv, h, x, S->parent, c);
        if (top < 0)
        {
            ok = top;
            break;
        }
        if (mpz_sgn(x[k]) != 0)
        {
            pivot = k;
            mpz_set(rhos[k], x[k]);
        }
        else
            pivot = SLIP_SINGULAR;
            
        // Error
        if (pivot < SLIP_OK)
        {
            ok = pivot;
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
                ok = slip_mpz_init2(L->x[lnz], size+2);
                if (ok != SLIP_OK) break;
                // Place the x value of the L->nz nonzero
                ok = slip_mpz_set(L->x[lnz],x[jnew]);
                if (ok != SLIP_OK) break;
                // Increment L->nz
                lnz += 1;
            }
        }
    }
    L->nz = S->lnz;

    //--------------------------------------------------------------------------
    // find the elasped time
    //--------------------------------------------------------------------------
    std::chrono::steady_clock::time_point t_end 
        = std::chrono::steady_clock::now();
    option->t_f = std::chrono::duration_cast<std::chrono::duration<float>>
        (t_end - t_begin);

    //--------------------------------------------------------------------------
    // Set the determinant of A (may be scaled)
    //--------------------------------------------------------------------------
    if (ok == SLIP_OK)
        ok = slip_mpz_set(option->determinant, rhos[n-1]);
    // Print statistics if desired
    if (option->print == 1)
    {
        std::cout<<"\n\n****Factorization Statistics****";
        std::cout<<"\n\nNum Nonzeros in A: "<<A->nz;
        std::cout<<"\nDimension of A, L, U: "<<n;
        std::cout<<"\nEstimated L nonzeros: "<<S->lnz;
        std::cout<<"\nActual L nonzeros: "<<L->nz;
        std::cout<<"\nNum Costly Reallocs: "<<numRealloc;
        std::cout<<"\n\n";
    }

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------
    SLIP_delete_mpz_array(x,n);
    SLIP_free(xi); SLIP_free(h); 
    // Collapse L
    //if (ok == SLIP_OK)
    //    SLIP_mat_collapse(L); */
    return ok;
}

#endif