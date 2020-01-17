#ifndef SLIP_Chol
#define SLIP_Chol

#include "../SLIP_LU/v0.6/Headers/SLIP_LU_config.h"
#include "UP_chol.h"



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

int LChol_ereach 
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
        if (i > k)
        {
            continue;
        }
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

/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 * This version uses left-looking LU's symbolic analysis
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
    top = cs_reach(L, A, k, xi, pinv);     // Obtain nonzero pattern in xi[top..n]
    
    //std::cout<<"\nReal Top is: "<<top;
    //std::cout<<"\nn-top is: "<<n-top;
    std::sort(xi + top, xi + n);
    //top = cs_reach(L, A, k, xi, pinv);     // Obtain nonzero pattern in xi[top..n]
    //SLIP_sort_xi(xi, top, n, pinv, pinv);  // Sort xi wrt sequence of pivots
    //SLIP_sort_xi(xi, top, n, pinv, row_perm);  // Sort xi wrt sequence of pivots
    SLIP_reset_mpz_array_2(x,n,top,xi);      // Reset x[i] = 0 for all i in nonzero pattern
    mpz_set_ui(x[k], 0);                    // x[k] may not be zeroed out properly!
    SLIP_reset_int_array_2(h,n,top,xi);      // Reset h[i] = -1 for all i in nonzero pattern
    SLIP_get_column(A,k,x);               // Set x = A(:,k)
    
    //--------------------------------------------------------------------------
    // Set the first k terms of x from row k of L
    //--------------------------------------------------------------------------
    
     for (j = top; j < n; j++)
     {
         p = xi[j];
         if (p < k)
         {
             for (i = L->p[p]; i < L->p[p+1]; i++)
             {
                 if (L->i[i] == k)
                 {
                     mpz_set(x[p], L->x[i]);
                 }
             }
         }
     }
    
    /*for (j = 0; j < k; j++)
    {
        for (i = L->p[j]; i < L->p[j+1]; i++)
        {
            if (L->i[i] == k)
            {
                mpz_set(x[j], L->x[i]);
                break;
            }
        }
    }*/

    //--------------------------------------------------------------------------
    // Iterate accross nonzeros in x
    //--------------------------------------------------------------------------
    ok = SLIP_OK;
    for (p = top; p < n; p++)
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
    //S->lnz = 3;
    //std::cout<<"\nS->lnz is: "<<S->lnz;
    ok = SLIP_mat_alloc2(n, n, S->lnz, L);
    //for (k = 0; k < n; k++) L->p[k] = c[k] = S->cp[k];
    //L->p[n] = S->lnz;
    if (ok != SLIP_OK) return ok;
    
    //--------------------------------------------------------------------------
    // Iteration 0, must select pivot
    //--------------------------------------------------------------------------
    check = 0;
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
        j+=1;
    }
    
    // Some error
    if (pivot < SLIP_OK)
        return pivot;
    std::sort(xi + top, xi + n);
    //SLIP_sort_xi(xi, top, n, pinv, pinv);  // Sort xi wrt sequence of pivots
    // Populate L and U
    for (j = top; j < n; j++)
    {
        jnew = xi[j];
        if (jnew >= 0)
        {
            // ith value of x[j]
            L->i[lnz] = jnew;
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
        // Column pointers for column k of L and U
        L->nz = lnz; 
        L->p[k] = lnz;
                       

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
            //std::cout<<"\n jnew is: "<<jnew;
            //std::cout<<"\nlnz is: "<<lnz;
            //std::cout<<"\nok is: "<<ok;
            //std::cout<<"\nk is: "<<k;
            if (jnew >= k && ok == SLIP_OK)
            {
                // Place the i location of the L->nz nonzero
                //std::cout<<"\nHey";
                L->i[lnz] = jnew;
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
    std::cout<<"\nAt end, lnz is: "<<lnz;
    L->nz = S->lnz;
    // Finalize L->p
    L->p[n] = S->lnz;
        

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


/* Preallocate the matrix L */


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
int SLIP_chol_triangular_solve2 // performs the sparse REF triangular solve
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
   
   // Assign nonzero pattern
     top = n - (L->p[k+1] - L->p[k]);
     int count = top;
     for (j = 0; j < k; j++)
     {
         for (i = L->p[j]; i < L->p[j+1]; i++)
         {
             if (L->i[i] == k)
             {
                 top-=1;
                 xi[top] = j;
                 mpz_set(x[j], L->x[i]);
                 break;
             }
         }
     }
     
     for (i = L->p[k]; i < L->p[k+1]; i++)
     {
         xi[count] = L->i[i];
         count+=1;
     }
    //top = Chol_ereach(A, k, parent, xi, c);
    std::sort(xi+top, xi+n);   
    SLIP_reset_mpz_array_2(x,n,top,xi);      // Reset x[i] = 0 for all i in nonzero pattern
    mpz_set_ui(x[k], 0);                    // x[k] may not be zeroed out properly!
    SLIP_reset_int_array_2(h,n,top,xi);      // Reset h[i] = -1 for all i in nonzero pattern
    SLIP_get_column(A,k,x);               // Set x = A(:,k)
    
     for (j = top; j < n; j++)
     {
         p = xi[j];
         if (p < k)
         {
             for (i = L->p[p]; i < L->p[p+1]; i++)
             {
                 if (L->i[i] == k)
                 {
                     mpz_set(x[p], L->x[i]);
                 }
             }
         }
     }
        
    
    //--------------------------------------------------------------------------
    // Iterate accross nonzeros in x
    //--------------------------------------------------------------------------
    ok = SLIP_OK;
    for (p = top; p < n; p++)
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
int SLIP_Chol_Factor2         // performs the SLIP LU factorization
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
        top = SLIP_chol_triangular_solve2(L, A, k, xi, rhos, pinv, h, x, S->parent, c);
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















/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
int SLIP_chol_triangular_solve3 // performs the sparse REF triangular solve
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
        //std::cout<<"\nc[m] is: "<<c[m];
        p = c[m]++;
        p+=1;
        //p++;
        //std::cout<<"\nnow c[m] is: "<<c[m];
        //std::cout<<"\nm is: "<<m;
       // std::cout<<"\np is: "<<p;
        mpz_set(x[m], L->x[p]);
        //std::cout<<"\nx["<<m<<"] is: "<<x[m];
        //mpz_set(x[m], L->x[p]);
        //std::cout<<"\np is: "<<p;
        //std::cout<<"\nL->x[p] is: "<<L->x[p];
        //std::cout<<"\nAt iteration "<<k<<" For nonzero "<<p << " the L index is: "<<L->i[ p] <<" or " <<L->i[p++];
    }

    //std::cout<<"\nAt iteration :" << k<< "\n";
    //for (i = top; i < n; i++)
    //    std::cout<<"\n x["<<xi[i]<<"] is: "<<x[xi[i]];
    //    //mpz_set( x [L->
    for (i = A->p[k]; i < A->p[k+1]; i++)
    {
        if ( A->i[i] >= k)
        {
            mpz_set(x[A->i[i]], A->x[i]);
        }
    }
    std::sort(xi+top, xi+n); 
    //if (k == 4)
    //{
        //std::cout<<"\nColumn 0 of L is: \n";
       // for (j = L->p[0]; j < L->p[1]; j++)
      //      std::cout<<"Entry in position: "<<L->i[j] <<" is : "<<L->x[j];
        
     //   std::cout<<"\nColumn 1 of L is: \n";
    //    for (j = L->p[1]; j < L->p[2]; j++)
    //        std::cout<<"Entry in position: "<<L->i[j] <<" is : "<<L->x[j];
    //}
    //if (k == 4) return -1;
   
    // FIGURE OUT HOW TO IMPLEMENT THIS!!!!!
    
    
    
    
    //SLIP_reset_mpz_array(x,n);      // Reset x[i] = 0 for all i in nonzero pattern
    //int* marked = SLIP_initialize_int_array(n);
    
    
    // xi is all nonzeros in location 0..k. I already know 
    
//     std::cout<<"\nAt iteration "<<k<< " the indices are: \n";
//     for (i = top; i < n; i++)
//         std::cout<< " " << xi[i];
//     std::cout<<"\nThe other nonzeros in this column are: \n";
//     for (i = L->p[k]; i < L->p[k+1]; i++)
//         std::cout<< " " << L->i[i];
    
    //std::cout<<"\nn is: "<<n;
    //std::cout<<"\nTop is: "<<top;
    //top--;
    //xi[top] = k;
    //for (i = 0; i < k; i++)
    //{
      //  for (j = L->p[i]; j < L->p[i+1]; j++)
//        {
  //          if (L->i[j] == k)
    //        {
      //          top--;
        ///        xi[top] = i;
           //     mpz_set(x[i], L->x[j]);
            //}
        //}
    //}
//     j = top;
//     for ( ; j < n; j++)
//     {
//         i = xi[j];
//         if (i > k) continue;
//         //if (i == k)
//         //    std::cout<<"\nHello";
//         marked[i] = 1;
//         p = c[i]++;
//         mpz_set(x[i], L->x[p++]);
//         for ( ; p < L->p[i+1]-1; p++)
//         {
//             if (L->i[p] < k && marked[i] < 1)
//             {
//                 top--;
//                 xi[top] = L->i[p];
//                 mpz_set(x[L->i[p]], L->x[p]);
//                 marked[L->i[p]] = 1;
//             }
//         }
//     }
    
    //for (p = A->p[k]; p < A->p[k+1]; p++)
    //{
      //  if (A->i[p] >= k)
//        {
  //          mpz_set(x[A->i[p]], A->x[p]);
    //    }
    //}
    //std::cout<<"\nHeyheyhey";
    //for (j = top; j < n; j++)
    //{
      //  if (xi[j] == k)
        //    std::cout<<"\nHere";
    //}
    //std::sort(xi+top, xi+n);   
    //std::cout<<"\nNonzero pattern is:\n";
    //for (j = top; j < n; j++)
    //    std::cout<<" "<<xi[j];
    //SLIP_reset_mpz_array_2(x,n,top,xi);      // Reset x[i] = 0 for all i in nonzero pattern
    //mpz_set_ui(x[k], 0);                    // x[k] may not be zeroed out properly!
    SLIP_reset_int_array_2(h,n,top,xi);      // Reset h[i] = -1 for all i in nonzero pattern
    //std::cout<<"\nx[k] is: "<<x[k];
    
    
    
    //SLIP_get_column(A,k,x);               // Set x = A(:,k)
    
//     for (j = top; j < n; j++)
//     {
//        i = xi[j];
//         if (i < k)
//         {
//            p = c[i]++;
//             for ( ; p < L->p[i+1]-1; p++)
//             {
//                mpz_set(x[L->i[p]], L->x[p]);
//             }
//         }
//     }
    
    /* for (j = top; j < n; j++)
     {
         p = xi[j];
         if (p < k)
         {
             for (i = L->p[p]; i < L->p[p+1]; i++)
             {
                 if (L->i[i] == k)
                 {
                     mpz_set(x[p], L->x[i]);
                 }
             }
         }
     }
    */    
    
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
int SLIP_Chol_Factor3         // performs the SLIP LU factorization
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
        top = SLIP_chol_triangular_solve3(L, A, k, xi, rhos, pinv, h, x, S->parent, c);
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