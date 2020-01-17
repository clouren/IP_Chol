//------------------------------------------------------------------------------
// IP_Chol/IP_tripread_double: Read in a triplet matrix
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* Purpose: This function reads in a double matrix stored in a triplet format
 * This format used can be seen in any of the example mat files. 
 * 
 * This is only used for Demo purposes
 */

SLIP_info IP_tripread_double
(
    SLIP_sparse* A,        // Matrix to be populated
    FILE* file,          // file to read from (must already be open)
    int32_t** i,
    int32_t** j,
    double** x,
    int* n,
    int* nz,
    SLIP_options* option
)
{
    int32_t m;
    SLIP_info ok;
    if (A == NULL || file == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }
    // Read in triplet form first
    
    // Read in size of matrix & number of nonzeros
    ok = fscanf(file, "%d %d %d\n", &m, n, nz);
    
    if (feof(file) || ok < 3)
    {
        return SLIP_INCORRECT_INPUT;
    }
    
    (*i) = (int32_t*) SLIP_malloc( (*nz) * sizeof(int32_t));
    (*j) = (int32_t*) SLIP_malloc( (*nz)* sizeof(int32_t));
    (*x) = (double*) SLIP_malloc( (*nz) * sizeof(double));

    if (!i || !j || !x )
    {
        return SLIP_OUT_OF_MEMORY;
    }

    int32_t decrement;
    ok = fscanf(file, "%d %d %lf\n", &( (*i) [0]), &( (*j) [0]), &( (*x) [0]));
    if (feof(file) || ok < 3)
    {
        return SLIP_INCORRECT_INPUT;
    }

    if (SLIP_MIN( (*i) [0], (*j) [0]) == 0)
    {
        decrement = 0;
    }
    else
    {
        decrement = 1;
        (*i)[0]-=decrement;
        (*j)[0]-=decrement;
    }

    // Read in the values from file
    for (int32_t k = 1; k < *nz; k++)
    {
        ok = fscanf(file, "%d %d %lf\n", &( (*i) [k]), &( (*j) [k]), &( (*x) [k]));
        if ((feof(file) && k != *nz-1) || ok < 3)
        {
            return SLIP_INCORRECT_INPUT;
        }
        // Conversion from 1 based to 0 based
        (*i)[k] -= decrement;
        (*j)[k] -= decrement;
    }
    return SLIP_OK;
}