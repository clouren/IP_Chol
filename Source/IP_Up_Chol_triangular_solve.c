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
int IP_Up_Chol_triangular_solve // performs the sparse REF triangular solve
(
    SLIP_sparse* L,              // partial L matrix
    SLIP_sparse* A,              // input matrix
    int k,                    // iteration of algorithm
    int* xi,                  // nonzero pattern vector
    int* parent,              // Elimination tree
    int* c,                   // Column pointers
    mpz_t* rhos,              // sequence of pivots
    int* h,                   // history vector
    mpz_t* x                  // solution of system ==> kth column of L and U
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
    IP_reset_mpz_array(x,n,top,xi);      // Reset x[i] = 0 for all i in nonzero pattern
    OK(SLIP_mpz_set_ui(x[k], 0));                     // Reset x[k]
    IP_reset_int_array2(h,n,top,xi);      // Reset h[i] = -1 for all i in nonzero pattern
    OK(IP_Up_get_column(A,k,x));                    // Set x = A(:,k)
    
    //--------------------------------------------------------------------------
    // Iterate accross nonzeros in x
    //--------------------------------------------------------------------------
    for (p = top; p < n; p++)
    {   
        /* Finalize x[j] */
        j = xi[p];                          // First nonzero term
        //if (j == k) continue;             // Do not operate on x[k] in TS
        if (mpz_sgn(x[j]) == 0) continue;   // If x[j] == 0 no work must be done
                
        // History update x[j]
        if (h[j] < j-1)
        {
            // History update x[j]
            OK(SLIP_mpz_mul(x[j], x[j], rhos[j-1]));
            if (h[j] > -1)
            {
               OK(SLIP_mpz_divexact(x[j], x[j], rhos[ h[j]]));
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
                        OK(SLIP_mpz_submul(x[i],L->x[m],x[j]));// x[i] = 0 - lij*x[j]
                        h[i] = j;                  // Entry is up to date
                    }
                        
                    // Previous pivot exists
                    else
                    {
                        OK(SLIP_mpz_submul(x[i],L->x[m],x[j]));// x[i] = 0 - lij*x[j]
                        OK(SLIP_mpz_divexact(x[i],x[i],rhos[j-1]));// x[i] = x[i] / rho[j-1]
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
                        OK(SLIP_mpz_mul(x[i],x[i],rhos[0]));      // x[i] = x[i]*rho[0]
                        OK(SLIP_mpz_submul(x[i], L->x[m], x[j]));// x[i] = x[i] - lij*xj
                        h[i] = j;                  // Entry is now up to date
                    }
                    // There is a previous pivot
                    else
                    {
                        // History update if necessary
                        if (h[i] < j - 1)
                        {
                            OK(SLIP_mpz_mul(x[i],x[i],rhos[j-1]));// x[i] = x[i] * rho[j-1]
                            if (h[i] > -1)
                            {
                                OK(SLIP_mpz_divexact(x[i],x[i],rhos[h[i]]));// x[i] = x[i] / rho[h[i]]
                            }
                        }
                        OK(SLIP_mpz_mul(x[i],x[i],rhos[j]));// x[i] = x[i] * rho[j]
                        OK(SLIP_mpz_submul(x[i], L->x[m], x[j]));// x[i] = x[i] - lij*xj
                        OK(SLIP_mpz_divexact(x[i],x[i],rhos[j-1]));// x[i] = x[i] / rho[j-1] 
                        h[i] = j;                  // Entry is up to date
                    }
                }
            }
        }
        // Update x[k]
        if (h[k] < j - 1)
        {
            OK(SLIP_mpz_mul(x[k],x[k],rhos[j-1])); // x[k] = x[k] * rho[j-1]
            if (h[k] > -1)
            {
                OK(SLIP_mpz_divexact(x[k],x[k],rhos[h[k]]));// x[k] = x[k] / rho[h[k]]
            }
        }
        OK(SLIP_mpz_mul(x[k],x[k],rhos[j])); // x[k] = x[k] * rho[j]
        OK(SLIP_mpz_submul(x[k], x[j], x[j]));// x[k] = x[k] - xj*xj
        if (j != 0)
            OK(SLIP_mpz_divexact(x[k],x[k],rhos[j-1])); // x[k] = x[k] / rho[j-1] 
        h[k] = j;   
    }
    // Finalize x[k]
    if (h[k] < k-1)
    {
        OK(SLIP_mpz_mul(x[k], x[k], rhos[k-1]));
        if (h[k] > -1)
        {
            OK(SLIP_mpz_divexact(x[k], x[k], rhos[ h[k]]));
        }
    }
    return top;
}