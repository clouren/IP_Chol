//------------------------------------------------------------------------------
// REF_Chol/IP_tripread_double: Read in a triplet matrix
//------------------------------------------------------------------------------

// REF Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// Texas A&M University.  All Rights Reserved.  See REF_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/REF-Chol.h"

/* Purpose: This function reads in a matrix stored in a triplet format
 * with double entries. The format used can be seen in any of the
 * example mat files.
 * 
 * This is only used for Demo purposes
 */

IP_Chol_info IP_tripread_double
(
    SLIP_matrix **A_handle,     // Matrix to be populated
    FILE* file,                 // file to read from (must already be open)
    SLIP_options* option        // Command options
)
{
    IP_Chol_info info ;
    if (A_handle == NULL || file == NULL)
    {
        printf ("invalid input\n") ;
        return IP_Chol_INCORRECT_INPUT;
    }
    (*A_handle) = NULL ;

    // Read in triplet form first
    int64_t m, n, nz;

    // Read in size of matrix & number of nonzeros
    int s = fscanf(file, "%"PRId64" %"PRId64" %"PRId64"\n", &m, &n, &nz);
    if (feof(file) || s < 3)
    {
        printf ("premature end-of-file\n") ;
        return SLIP_INCORRECT_INPUT;
    }

    // First, we create our A matrix which is triplet double
    SLIP_matrix *A = NULL;
    info = SLIP_matrix_allocate(&A, SLIP_TRIPLET, SLIP_FP64, m, n, nz,
        false, true, option);
    if (info != IP_Chol_OK)
    {
        return (info) ;
    }
    
    s = fscanf (file, "%"PRId64" %"PRId64" %lf\n",
        &(A->i[0]), &(A->j[0]), &(A->x.fp64[0])) ;
            
    if (feof(file) || s <= 0)
    {
        printf ("premature end-of-file\n") ;
        SLIP_matrix_free(&A, option);
        return SLIP_INCORRECT_INPUT;
    }

    // Matrices in this format are 1 based. We decrement
    // the indices by 1 to use internally
    A->i[0] -= 1;
    A->j[0] -= 1;

    // Read in the values from file
    for (int64_t k = 1; k < nz; k++)
    {
        s = fscanf(file, "%"PRId64" %"PRId64" %lf\n",
            &(A->i[k]), &(A->j[k]), &(A->x.fp64[k]));
        if ((feof(file) && k != nz-1) || s < 3)
        {
            printf ("premature end-of-file\n") ;
            SLIP_matrix_free(&A, option);
            return SLIP_INCORRECT_INPUT;
        }
        // Conversion from 1 based to 0 based
        A->i[k] -= 1;
        A->j[k] -= 1;
    }

    // the triplet matrix now has nz entries
    A->nz = nz;    

    // At this point, A is a double triplet matrix. We make a copy of it with C
    // C is a CSC matrix with mpz entries
    //option->print_level = 3;
   // SLIP_matrix_check(A,option);
    
    SLIP_matrix* C = NULL;
    SLIP_matrix_copy(&C, SLIP_CSC, SLIP_MPZ, A, option);

    // Success. Set A_handle = C and free A

    SLIP_matrix_free(&A, option);
    (*A_handle) = C;
    return (info) ;
}
