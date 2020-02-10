///------------------------------------------------------------------------------
// IP_Chol/IP_determine_symmetry: This function determines if the input matrix is symmetric.
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Determine if the input A is indeed symmetric prior to factorization.
 * There are two options as to how to determine the symmetry. 
 * By setting the input exhaustive = 1, both the nonzero pattern and the values
 * of the nonzero entries are checked for symmetry. If A passes both of these tests,
 * then we can be sure it is indeed fully symmetric.
 * 
 * If exhaustive is set to any other value, only the nonzero pattern of A is checked,
 * thus we cannot gauranteee that the matrix is indeed fully symmetric as the values
 * of the entries is not checked.
 * 
 * On success, 0 is returned. If the matrix is not symmetric, 1 is returned.
 * 
 */

#include "../Include/IP-Chol.h"

int IP_determine_symmetry
(
    SLIP_sparse* A,
    int exhaustive
)
{
    int j;
    SLIP_sparse* T = IP_transpose(A);
    for (j = 0; j < A->nz; j++)
    {
        if (T->i[j] != A->i[j])
        {
            printf("\nError, matrix is not symmetric\n");
            SLIP_delete_sparse(&T);
            return 1;
        }
    }
    
    for (j = 0; j <= A->n; j++)
    {
        if (T->p[j] != A->p[j])
        {
            printf("\nError, matrix is not symmetric\n");
            SLIP_delete_sparse(&T);
            return 1;
        }
    }
    
    if (exhaustive == 1)
    {
        int r;
        for (j = 0; j < A->nz; j++)
        {
            SLIP_mpz_cmp(&r, A->x[j], T->x[j]);
            if ( r != 0)
            {
                printf("\nError, matrix is not symmetric\n");
                SLIP_delete_sparse(&T);
                return 1;
            }
        }
    }
    SLIP_delete_sparse(&T);
    return 0;
        
}