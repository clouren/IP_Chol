//------------------------------------------------------------------------------
// IP_Chol/IP_Up_Chol_triangular_solve: Sparse symmetric REF Triangular solve
//------------------------------------------------------------------------------

// IP Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// and Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"


/* Purpose: This function performs the symmetric sparse REF triangular solve for
 * the up-looking Cholesky factorization. i,e, LD x = A(1:k-1, k). At the end of
 * this function, the vector x contains the values of the kth row of the intger-
 * preserving matrix L. 
 * 
 * Command input:
 * top_output:      A pointer to the beginning of the nonzero pattern. Undefined
 *                  on input, on output xi[top_output..n] contains the beginning
 *                  of the nonzero pattern.
 * 
 * L:               Lower triangular matrix.
 * 
 * A:               Input matrix
 * 
 * k:               Current iteration of the algorithm
 * 
 * xi:              Nonzero pattern. Undefined on input, on output contains teh 
 *                  nonzero pattern of the kth row of L
 * 
 * parent:          Elimintaation tree
 * 
 * c:               Column pointers of L
 * 
 * rhos:            Pivot matrix
 * 
 * h:               History vector
 * 
 * x:               Solution of linear system. Undefined on input, on output
 *                  contains the kth row of L.
 */

// TODO Properly set which of these input/output are const

IP_Chol_info IP_Up_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,            // Output the beginng of nonzero pattern
    SLIP_matrix* L,                 // partial L matrix
    SLIP_matrix* A,                 // input matrix
    int64_t k,                      // iteration of algorithm
    int64_t*  xi,                   // nonzero pattern vector
    int64_t*  parent,               // Elimintaation tree
    int64_t* c,                     // Column point64_t*ers
    SLIP_matrix* rhos,              // sequence of pivots
    int64_t * h,                    // history vector
    SLIP_matrix* x                  // solution of system ==> kth column of L and U
)
{
    IP_Chol_info ok;
    
    SLIP_REQUIRE(L, SLIP_CSC, SLIP_MPZ);
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);
    SLIP_REQUIRE(rhos, SLIP_DENSE, SLIP_MPZ);
    SLIP_REQUIRE(x, SLIP_DENSE, SLIP_MPZ);
    
    if (!xi || !parent || !c || !h)
        return SLIP_INCORRECT_INPUT;
    
    int64_t j, i, inew, p, m, top, n = A->n, col;
    
    //--------------------------------------------------------------------------
    // Initialize REF TS by getting nonzero patern of x && obtaining A(:,k)
    //--------------------------------------------------------------------------
    
    // Obtain the nonzero pattern of row k of L. The nonzero pattern is stored in
    // xi[top..n]
    top = IP_Chol_ereach(A, k, parent, xi, c);
    
    // Sort the nonzero pattern of row k of L
    qsort(&xi[top], n-top, sizeof(int64_t*), compare); 
    
    
    // Reset x[i] = 0 for all i in nonzero pattern xi [top..n-1]
    for (i = top; i < n; i++)
    {
        SLIP_CHECK (SLIP_mpz_set_ui (x->x.mpz[xi [i]], 0)) ;
    }
    
    // Reset value of x[k]. If the matrix is nonsingular, x[k] will
    // be a part of the nonzero pattern and reset in the above loop.
    // If the matrix is detected to be singular at this iteration,
    // x[k] will be 0. However, in some rare cases, the matrix can be singular 
    // but x[k] will be nonzero from a previous iteration (since entries are 
    // only overwritten when necessary). Thus, here we manually reset
    // x[k] to account for this extremely rare case.
    SLIP_CHECK( SLIP_mpz_set_ui( x->x.mpz[k], 0));
        
    // Reset the history vector, h[i] = -1 for all i in nonzero pattern
    for (i = top; i < n; i++)
    {
        h[xi[i]] = -1;
    }
    
    // Initialization, Set x = A(:,k)
    for (i = A->p[k]; i < A->p[k+1]; i++)
    {
        // Only use the nonzeros in the upper triangular portion of 
        // column k of A
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
        j = xi[p];                                  // First nonzero term

        if (mpz_sgn(x->x.mpz[j]) == 0) continue;    // If x[j] == 0 no work must be done
                
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
    
    // At this point in the algorithm every entry in L(k,:) is finalized except for potentially 
    // the pivot element. It is possible that the pivot element is not at its final value,
    // thus we perform a final history update to finalize entry x[k]
    
    if (h[k] < k-1)
    {
        OK(SLIP_mpz_mul(x->x.mpz[k], x->x.mpz[k], rhos->x.mpz[k-1]));
        if (h[k] > -1)
        {
            OK(SLIP_mpz_divexact(x->x.mpz[k], x->x.mpz[k], rhos->x.mpz[ h[k]]));
        }
    }
    // Set the beginning of the nonzero pattern.
    *top_output = top;
    return SLIP_OK;
}
