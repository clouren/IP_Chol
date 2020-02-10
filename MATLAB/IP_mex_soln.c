//------------------------------------------------------------------------------
// IP-Chol/MATLAB/IP_Chol_mex_soln: Use IP-Chol within MATLAB
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "./Source/IP-Chol_mex.h"
//
//#include "../SLIP_LU-master/SLIP_LU/Include/SLIP_LU.h"


/* Purpose: .c files defining the IP-Chol matlab interfacee
 * This function defines: x = IP_Chol(A, b, option)
 */


void mexFunction
(
    int32_t nargout,
    mxArray *pargout [ ],
    int32_t nargin,
    const mxArray *pargin [ ]
)
{
    //--------------------------------------------------------------------------
    // Initialize SLIP LU library environment
    //--------------------------------------------------------------------------
    SLIP_initialize_expert(mxMalloc, IP_gmp_mex_realloc, IP_gmp_mex_free);
    SLIP_info status;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    IP_check_input(pargin, nargin);
    if (nargout > 1 || nargout <= 0 || nargin != 3)
    {
        mexErrMsgTxt("Usage: x = SLIP_LU(A,b,option)");
    }
    //--------------------------------------------------------------------------
    // Allocate memory
    //--------------------------------------------------------------------------
    SLIP_sparse *A = NULL, *L = NULL;
    A = SLIP_create_sparse();
    L = SLIP_create_sparse();
    SLIP_dense *b = SLIP_create_dense();
    //Set defaults for options
    SLIP_options* option = SLIP_create_default_options();
    if (!A || !L || !b || !option)
    {
        IP_mex_error (SLIP_OUT_OF_MEMORY);
    }
    SLIP_LU_analysis* S = NULL;

    //--------------------------------------------------------------------------
    // Declare variables and process input
    //--------------------------------------------------------------------------
    // Read in options
    IP_get_matlab_options(option, pargin[2]);

    // Read in A and b
    IP_mex_get_A_and_b(A, b, pargin, nargin);

    // Create arrays based on the size of input matrix
    S = SLIP_create_LU_analysis((A->n)+1);
    double** soln = SLIP_create_double_mat(b->m, b->n);
    int32_t* pinv = (int32_t*) SLIP_malloc(A->n* sizeof(int32_t));
    mpz_t* rhos = SLIP_create_mpz_array(A->n);
    mpq_t** soln_mpq = SLIP_create_mpq_mat(b->m, b->n);
    if (!S || !soln || !pinv || !rhos || !soln_mpq)
    {
        IP_mex_error (SLIP_OUT_OF_MEMORY);
    }

    option->order = SLIP_AMD;  // AMD
    //--------------------------------------------------------------------------
    // Symbolic analysis and factorization
    //--------------------------------------------------------------------------
    
    
    IP_MEX_OK (SLIP_LU_analyze(S, A, option));// Symbolic Analysis

    int* pinv2 = (int*) SLIP_malloc(A->n* sizeof(int));
    for (int k = 0; k < A->n; k++)
    {
        int index = S->q[k];
        pinv2[index] = k;
    }
    
    SLIP_sparse* A2 = NULL;
    A2 = IP_Chol_permute_A(A, pinv2, S);
    Sym_chol* S2 = (Sym_chol*) SLIP_malloc(1* sizeof(Sym_chol));
    IP_MEX_OK(IP_Up_Chol_Factor(A2, L, S2, rhos, option));
    
    //--------------------------------------------------------------------------
    // FB Substitution
    //--------------------------------------------------------------------------
    
    IP_MEX_OK(IP_Solve(soln_mpq, b->x, rhos, L, pinv2, option, 1));
        
    IP_MEX_OK(SLIP_permute_x(soln_mpq, b->m, b->n, S));

    IP_MEX_OK(SLIP_scale_x(soln_mpq, A, b));
    IP_MEX_OK(SLIP_get_double_soln(soln, soln_mpq, b->m, b->n));

    //--------------------------------------------------------------------------
    // Set outputs, free memory
    //--------------------------------------------------------------------------
    pargout[0] =  IP_mex_output_soln(soln, b->m, b->n);
    SLIP_delete_mpq_mat(&soln_mpq, b->m, b->n);
    SLIP_delete_mpz_array(&rhos, A->n);
    SLIP_FREE(pinv);
    SLIP_delete_double_mat(&soln, b->m, b->n);
    SLIP_delete_LU_analysis(&S);
    SLIP_FREE(option);
    SLIP_delete_dense(&b);
    SLIP_delete_sparse(&A2);
    SLIP_delete_sparse(&L);
    SLIP_delete_sparse(&A);
    SLIP_finalize();
}
