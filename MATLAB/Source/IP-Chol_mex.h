//------------------------------------------------------------------------------
// IP-Chol/MATLAB/IP-Chol_mex.h: Include file for IP-Chol
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#ifndef IP_mex 
#define IP_mex

#include "../../Include/IP-Chol.h"
# include "mex.h"
#include "matrix.h"


#define IP_MEX_OK(method)         \
{                                   \
    status = method;                \
    IP_mex_error(status);         \
}


/* Purpose: A GMP reallocation function 
 * This allows GMP to use MATLAB's default realloc function 
 */
void* IP_gmp_mex_realloc 
(
    void* x,    // void* to be reallocated 
    size_t a,   // Previous size
    size_t b    // New size
);

/* Purpose: A GMP free function. This allows GMP to use
 * MATLAB's mxFree instead of free 
 */
void IP_gmp_mex_free 
(
    void* x,    // void* to be freed
    size_t a    // Size
);

/* Purpose: This function converts mpq array to double
 * NOTE: This induces roundoff error via the final division
*/
void IP_mpq_to_double
(
    double* x_doub,       // double array
    const mpq_t* x_mpq,   // mpq array
    const int32_t n       // size of b
);

void IP_check_input
(
    const mxArray * input [],    // The matlab input array
    int32_t nargin
);

void IP_get_matlab_options
(
    SLIP_options* option,  // Control parameters
    const mxArray* input   // The input options from MATLAB interface
);

/* Purpose: Convert int32_t* array to MATLAB mwIndex* array
 */
void IP_int32_to_mwIndex
(
    mwIndex* y, 
    int32_t* x, 
    int32_t n
) ;

/* Purpose: convert mwIndex* array to int32_t* array
 */
void IP_mwIndex_to_int32
(
    int32_t* y, 
    mwIndex* x, 
    mwSize n
) ;

/* Purpose: Check the input x array for numbers too large for 
 * double precision.
 */
void IP_mex_check_for_inf
(
    double* x, // The array of numeric values 
    mwSize n   // size of array
);

/* Purpose: This function reads in the A matrix and right hand side vectors. */
void IP_mex_get_A_and_b
(
    SLIP_sparse *A,          // Internal SLIP Mat stored in ccf 
    SLIP_dense *b,           // mpz matrix used internally
    const mxArray* input[],  // The input A matrix and options 
    int32_t nargin           // Number of input to the mexFunction
);


/* Purpose: Output the solution to the linear system Ax=b to matlab
 */
mxArray* IP_mex_output_soln
(
    double** x, 
    int32_t m, 
    int32_t n
) ;


/* Purpose: Report errors if they arise
 */
void IP_mex_error
(
    SLIP_info status
) ;

/* Purpose: Drop entries which are zero from a sparse matrix
 */
mwIndex IP_dropzeros 
(
    mxArray *A
);

/* Purpose: Used to drop zeros 
 */
mwIndex IP_fkeep 
(
    mxArray *A, 
    bool (*fkeep) (int32_t, int32_t, double)
);

/* Purpose: transpose the matrix A
 */
SLIP_info IP_transpose_mex 
(
    mxArray *A
);

#endif
