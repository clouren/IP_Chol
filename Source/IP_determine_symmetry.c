///------------------------------------------------------------------------------
// IP_Chol/IP_determine_symmetry: This function determines if the input matrix is symmetric.
//------------------------------------------------------------------------------

// IP Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// and Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// IP_Chol/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Determine if the input A is indeed symmetric prior to factorization.
 * There are two options as to how to determine the symmetry. 
 * By setting the input exhaustive to true, both the nonzero pattern and the values
 * of the nonzero entries are checked for symmetry. If A passes both of these tests,
 * then we can be sure it is indeed fully symmetric.
 * 
 * If exhaustive is set to false, only the nonzero pattern of A is checked,
 * thus we cannot gauranteee that the matrix is indeed fully symmetric as the values
 * of the entries is not checked.
 * 
 * On success, 0 is returned. If the matrix is not symmetric, 1 is returned.
 * 
 */

#include "../Include/IP-Chol.h"

int64_t IP_determine_symmetry
(
    SLIP_matrix* A,
    bool exhaustive
)
{
    int64_t j;
    
    // Declare matrix T
    SLIP_matrix *T = NULL;    
    // T = A'
    IP_transpose(&T, A);
    
    // Check if i values are the same
    for (j = 0; j < A->nz; j++)
    {
        if (T->i[j] != A->i[j])
        {
            printf("\nError, matrix is not symmetric\n");
            SLIP_matrix_free(&T,NULL);
            return 1;
        }
    }
    
    // Check if column pointers are the same
    for (j = 0; j <= A->n; j++)
    {
        if (T->p[j] != A->p[j])
        {
            printf("\nError, matrix is not symmetric\n");
            SLIP_matrix_free(&T,NULL);
            return 1;
        }
    }

    // If we are performing an exhaustive search, we check the x values as well
    // This is by far the most expensive part of checking the symmetry.
    if (exhaustive == true)
    {
        int r;
        for (j = 0; j < A->nz; j++)
        {
            SLIP_mpz_cmp(&r, A->x.mpz[j], T->x.mpz[j]);
            if ( r != 0)
            {
                printf("\nError, pattern is symmetric, values are not\n");
                SLIP_matrix_free(&T,NULL);
                return 1;
            }
        }
    }
    SLIP_matrix_free(&T,NULL);
    return 0;
        
}
