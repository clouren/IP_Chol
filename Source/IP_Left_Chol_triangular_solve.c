//------------------------------------------------------------------------------
// IP_Chol/IP_Left_Chol_triangular_solve: sparse symmetric left-looking triangular solve
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

static inline int compare4 (const void * a, const void * b)
{
    return ( *(int64_t*)a - *(int64_t*)b );
}


/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
SLIP_info IP_Left_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,        // Output the beginning of nonzero pattern
    SLIP_matrix* L,              // partial L matrix
    SLIP_matrix* A,              // input matrix
    int64_t k,                    // iteration of algorithm
    int64_t* xi,                  // nonzero pattern vector
    SLIP_matrix* rhos,              // sequence of pivots
    int64_t* h,                   // history vector
    SLIP_matrix* x,                  // solution of system ==> kth column of L and U
    int64_t* parent,
    int64_t* c
)
{
    SLIP_info ok;
    if (!L || !A || !xi || !rhos || !h || !x)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int64_t j, jnew, i, inew, p, m, top, n, col;

    //--------------------------------------------------------------------------
    // Initialize REF TS by getting nonzero patern of x && obtaining A(:,k)
    //--------------------------------------------------------------------------
    n = A->n;                                // Size of matrix and the dense vectors
    /* Chol_ereach gives the nonzeros located in L(k,:) upon completion
     * the vector xi contains the indices of the first k-1 nonzeros in column
     * k of L 
     */
    top = IP_Chol_ereach(A, k, parent, xi, c);
    j = top; // Store where the first k-1 nonzeros end
    // Populate the rest of the nonzero pattern
    for (i = L->p[k]; i < L->p[k+1]; i++)
    {
        top -= 1;           // One more nonzero in column k
        xi[top] = L->i[i];  // Index of the new nonzero
    }
    // Reset the array
    for (i = top; i < n; i++)
    {
        SLIP_mpz_set_ui(x->x.mpz[i], 0);
    }
    
    // Now we obtain the values of the first k-1 entries of x
    for (i = j; i < n; i++)
    {
        m = xi[i];
        p = c[m]++;
        p+=1;
        mpz_set(x->x.mpz[m], L->x.mpz[p]);
    }

    for (i = A->p[k]; i < A->p[k+1]; i++)
    {
        if ( A->i[i] >= k)
        {
            OK(SLIP_mpz_set(x->x.mpz[A->i[i]], A->x.mpz[i]));
        }
    }
    qsort(&xi[top], n-top, sizeof(int64_t), compare4);
    // Reset h[i] = -1 for all i in nonzero pattern
    for (i = top; i < n; i++)
    {
        h[xi[i]] = -1;
    }
        
    //--------------------------------------------------------------------------
    // Iterate accross nonzeros in x
    //--------------------------------------------------------------------------
    for ( p = top; p < n; p++)
    {   
        /* Finalize x[j] */
        j = xi[p];                        // First nonzero term
        if (mpz_sgn(x->x.mpz[j]) == 0) continue;// If x[j] == 0 no work must be done
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
                            OK(SLIP_mpz_submul(x->x.mpz[i], L->x.mpz[m], x->x.mpz[j]));// x[i] = x[i] - lij*xj
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
        }
        else                                              // Entries of L
        {
            //------------------------------------------------------------------
            // History update
            //------------------------------------------------------------------
            if (h[j] < k-1)
            {
                OK(SLIP_mpz_mul(x->x.mpz[j],x->x.mpz[j],rhos->x.mpz[k-1]));           // x[j] = x[j] * rho[k-1]
                if (h[j] > -1)
                {
                    OK(SLIP_mpz_divexact(x->x.mpz[j],x->x.mpz[j],rhos->x.mpz[h[j]]));// x[j] = x[j] / rho[h[j]]
                }
            }
        }
    }
    // Output the beginning of nonzero pattern
    *top_output = top;
    return SLIP_OK;
}
