//------------------------------------------------------------------------------
// IP_Chol/IP_Up_Chol_triangular_solve: Sparse symmetric REF Triangular solve
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"


static inline int32_t compare2 (const void * a, const void * b)
{
    return ( *(int32_t*)a - *(int32_t*)b );
}

/* Purpose: This function performs the symmetric sparse REF triangular solve. for uplooking
 * Cholesky factorization. i.e., 
 * (LD) x = A(1:k-1,k). 
 */
SLIP_info IP_Up_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,        // Output the beginning of nonzero pattern
    SLIP_matrix* L,              // partial L matrix
    SLIP_matrix* A,              // input matrix
    int k,                    // iteration of algorithm
    int* xi,                  // nonzero pattern vector
    int* parent,              // Elimination tree
    int* c,                   // Column pointers
    SLIP_matrix* rhos,              // sequence of pivots
    int* h,                   // history vector
    SLIP_matrix* x                  // solution of system ==> kth column of L and U
)
{
    SLIP_info ok;
    if (!L || !A || !xi || !parent || !c || !rhos || !h || !x)
        return SLIP_INCORRECT_INPUT;
    int j, jnew, i, inew, p, m, top, n = A->n, col;
    
    //--------------------------------------------------------------------------
    // Initialize REF TS by getting nonzero patern of x && obtaining A(:,k)
    //--------------------------------------------------------------------------
    top = IP_Chol_ereach(A, k, parent, xi, c);  // Obtain nonzero pattern in xi[top..n]
    qsort(&xi[top], n-top, sizeof(int), compare2); 
    
    
    // Reset x[i] = 0 for all i in nonzero pattern xi [top..n-1]
    for (i = top; i < n; i++)
    {
        SLIP_CHECK (SLIP_mpz_set_ui (x->x.mpz[xi [i]], 0)) ;
    }s
    OK(SLIP_mpz_set_ui(x[k], 0));                     // Reset x[k]
    
    
    // Reset h[i] = -1 for all i in nonzero pattern
    for (i = top; i < n; i++)
    {
        h[xi[i]] = -1;
    }
    
    // Set x = A(:,k)
    for (i = A->p[k]; i < A->p[k+1]; i++)
    {
        if (A->i[i] <= k)
        {
            OK(SLIP_mpz_set(x->x.mpz[A->i[i]], A->x.mpz[i]));
        }
    }
    
    //--------------------------------------------------------------------------
    // Iterate accross nonzeros in x
    //--------------------------------------------------------------------------
    for (p = top; p < n; p++)
    {   
        /* Finalize x[j] */
        j = xi[p];                          // First nonzero term
        //if (j == k) continue;             // Do not operate on x[k] in TS
        if (mpz_sgn(x->x.mpz[j]) == 0) continue;   // If x[j] == 0 no work must be done
                
        // History update x[j]
        if (h[j] < j-1)
        {
            // History update x[j]
            OK(SLIP_mpz_mul(x->x.mpz[j], x->x.mpz[j], rhos->x.mpz[j-1]));
            if (h[j] > -1)
            {
               OK(SLIP_mpz_divexact(x->x.mpz[j], x->x.mpz[j], rhos->x.mpz[ h[j]]));
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
                if (mpz_sgn(L->x.mpz[m]) == 0) continue;

                //----------------------------------------------------------
                /************* lij is nonzero, x[i] is zero****************/
                // x[i] = 0 then only perform IPGE update subtraction/division
                //----------------------------------------------------------
                if (mpz_sgn(x->x.mpz[i]) == 0)
                {
                    // No previous pivot
                    if (j < 1)
                    {
                        OK(SLIP_mpz_submul(x->x.mpz[i],L->x.mpz[m],x->x.mpz[j]));// x[i] = 0 - lij*x[j]
                        h[i] = j;                  // Entry is up to date
                    }
                        
                    // Previous pivot exists
                    else
                    {
                        OK(SLIP_mpz_submul(x->x.mpz[i],L->x.mpz[m],x->x.mpz[j]));// x[i] = 0 - lij*x[j]
                        OK(SLIP_mpz_divexact(x->x.mpz[i],x->x.mpz[i],rhos->x.mpz[j-1]));// x[i] = x[i] / rho[j-1]
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
                        OK(SLIP_mpz_mul(x->x.mpz[i],x->x.mpz[i],rhos->x.mpz[0]));      // x[i] = x[i]*rho[0]
                        OK(SLIP_mpz_submul(x->x.mpz[i], L->x->x.mpz[m], x->x.mpz[j]));// x[i] = x[i] - lij*xj
                        h[i] = j;                  // Entry is now up to date
                    }
                    // There is a previous pivot
                    else
                    {
                        // History update if necessary
                        if (h[i] < j - 1)
                        {
                            OK(SLIP_mpz_mul(x->x.mpz[i],x->x.mpz[i],rhos->x.mpz[j-1]));// x[i] = x[i] * rho[j-1]
                            if (h[i] > -1)
                            {
                                OK(SLIP_mpz_divexact(x->x.mpz[i],x->x.mpz[i],rhos->x.mpz[h[i]]));// x[i] = x[i] / rho[h[i]]
                            }
                        }
                        OK(SLIP_mpz_mul(x->x.mpz[i],x->x.mpz[i],rhos->x.mpz[j]));// x[i] = x[i] * rho[j]
                        OK(SLIP_mpz_submul(x->x.mpz[i], L->x.mpz[m], x->x.mpz[j]));// x[i] = x[i] - lij*xj
                        OK(SLIP_mpz_divexact(x->x.mpz[i],x->x.mpz[i],rhos->x.mpz[j-1]));// x[i] = x[i] / rho[j-1] 
                        h[i] = j;                  // Entry is up to date
                    }
                }
            }
        }
        // Update x[k]
        if (h[k] < j - 1)
        {
            OK(SLIP_mpz_mul(x->x.mpz[k],x->x.mpz[k],rhos->x.mpz[j-1])); // x[k] = x[k] * rho[j-1]
            if (h[k] > -1)
            {
                OK(SLIP_mpz_divexact(x->x.mpz[k],x->x.mpz[k],rhos->x.mpz[h[k]]));// x[k] = x[k] / rho[h[k]]
            }
        }
        OK(SLIP_mpz_mul(x->x.mpz[k],x->x.mpz[k],rhos->x.mpz[j])); // x[k] = x[k] * rho[j]
        OK(SLIP_mpz_submul(x->x.mpz[k], x->x.mpz[j], x->x.mpz[j]));// x[k] = x[k] - xj*xj
        if (j != 0)
            OK(SLIP_mpz_divexact(x->x.mpz[k],x->x.mpz[k],rhos->x.mpz[j-1])); // x[k] = x[k] / rho[j-1] 
        h[k] = j;   
    }
    // Finalize x[k]
    if (h[k] < k-1)
    {
        OK(SLIP_mpz_mul(x->x.mpz[k], x->x.mpz[k], rhos->x.mpz[k-1]));
        if (h[k] > -1)
        {
            OK(SLIP_mpz_divexact(x->x.mpz[k], x->x.mpz[k], rhos->x.mpz[ h[k]]));
        }
    }
    *top_output = top;
    return SLIP_OK;
}
