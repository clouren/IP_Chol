//------------------------------------------------------------------------------
// IP_Chol/SLIP_Chol: Solve an SPD linear system using left-chol
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

/* Include the Integer-preserving Cholesky routines */

#define FREE_WORKSPACE                  \
{                                       \
    SLIP_matrix_free(&A,NULL);          \
    SLIP_matrix_free(&L,NULL);          \
    SLIP_matrix_free(&A2,NULL);         \
    SLIP_FREE(pinv2);                   \
    SLIP_matrix_free(&b,NULL);          \
    SLIP_matrix_free(&rhos,NULL);       \
    SLIP_matrix_free(&x,NULL);          \
    SLIP_FREE(option);                  \
    SLIP_FREE(pinv);                    \
    SLIP_delete_LU_analysis(&S);        \
    SLIP_FREE(S2->parent);              \
    SLIP_FREE(S2->cp);                  \
    SLIP_FREE(S2);                      \
    SLIP_finalize();                    \
}                                       \

#define DEMO_OK(method)                 \
{                                       \
    ok = method ;                       \
    if (ok != SLIP_OK)                  \
    {                                   \
        IP_determine_error(ok);         \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}


# include "./Include/IP-Chol.h"


int main( int argc, char* argv[] )
{

    //--------------------------------------------------------------------------
    // Prior to using IP-Chol, its environment must be initialized. This is done
    // by calling the SLIP_initialize() function. 
    //--------------------------------------------------------------------------
    SLIP_initialize();
    
    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int n = 0, check, ok, j, index, k, nz = 0;
   
    SLIP_LU_analysis* S = NULL;
    SLIP_matrix *A = NULL;
    SLIP_matrix *L = NULL;
    SLIP_matrix *b = NULL;
    SLIP_matrix *rhos = NULL;
    int* pinv = NULL;
    int* pinv2 = NULL;
    SLIP_matrix* A2 = NULL;
    Sym_chol* S2 = NULL;
    SLIP_matrix* x = NULL;
    
    // Default options. May be changed in SLIP_LU_config.h
    SLIP_options *option = SLIP_create_default_options();
    
    char* mat_name = "./ExampleMats/2.mat";// Set demo matrix and RHS name
    char* rhs_name = "./ExampleMats/2.mat.soln";
    int rat = 1;
    
    // Process the command line
    DEMO_OK(IP_process_command_line(argc, argv, option,
        &mat_name, &rhs_name, &rat));
    
    //--------------------------------------------------------------------------
    // Allocate memory    
    //--------------------------------------------------------------------------
    
    // Read in A
    FILE* mat_file = fopen(mat_name,"r");
    if( mat_file == NULL )
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    
    // TODO Fix this function
    DEMO_OK(IP_tripread_double(A, mat_file, &I_in, &J_in, &x_in, &n, &nz, option));
    
    // For this code, we utilize a vector of all ones as the RHS vector    
    SLIP_matrix_allocate(&b, SLIP_DENSE, SLIP_MPZ, n, 1, n, false, true, option);
    pinv = (int*) SLIP_malloc(n* sizeof(int));
    // Create RHS
    for (int k = 0; k < n; k++)
        OK(SLIP_mpz_set_ui(b->x.mpz[k],1));
    
    //--------------------------------------------------------------------------
    // Perform Ordering of A
    //--------------------------------------------------------------------------
    clock_t start_col = clock();
        
    // Symmetric ordering of A. Uncomment the desired one, AMD is recommended
    //option->order = SLIP_NO_ORDERING;  // No ordering
    option->order = SLIP_AMD;  // AMD
    //option->order = SLIP_COLAMD; // COLAMD
        
    DEMO_OK(SLIP_LU_analyze(&S, A, option));
    
    clock_t end_col = clock();
    
    //--------------------------------------------------------------------------
    // Determine if A is indeed symmetric. If so, we try Cholesky
    // uncomment the one desired.
    // --------------------------------------------------------------------------
    
    clock_t start_sym = clock();
    
    int test = 0;
    //test = IP_determine_symmetry(A, 0);    // Determine symmetry just with nonzero pattern
    test = IP_determine_symmetry(A, 1);    // Determine symmetry with nonzero pattern and values
        
    if (test == 1) return 0;
    
    clock_t end_sym = clock();
    
    //--------------------------------------------------------------------------
    // Permute matrix A, that is set A2 = PAP'
    //--------------------------------------------------------------------------
    pinv2 = (int*) SLIP_malloc(n* sizeof(int));
    for (k = 0; k < n; k++)
    {
        index = S->q[k];
        pinv2[index] = k;
    }
    // TODO no return argument, fixme
    A2 = IP_Chol_permute_A(A, pinv2, S);
    
    //--------------------------------------------------------------------------
    // SLIP Chol Factorization
    //--------------------------------------------------------------------------
    clock_t start_factor = clock();
    
    // FIXME should be
    //OK(IP_Left_Chol_Factor(A2, &L, &S2, &rhos, &pinv, option));
    S2 = (Sym_chol*) SLIP_malloc(1* sizeof(Sym_chol));
    OK(IP_Left_Chol_Factor(A2, L, S2, rhos, pinv, option));
        
    L->m = n;
    clock_t end_factor = clock();
    
    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    clock_t start_solve = clock();
    // FIXME, mimic SLIP LU
    x = SLIP_create_mpq_mat(n, 1);
    DEMO_OK(IP_Solve(x, b->x, rhos, L, pinv2, option, 1));
    
    clock_t end_solve = clock();
    
    //--------------------------------------------------------------------------
    // Soln verification
    //--------------------------------------------------------------------------
    // x = Q x
    DEMO_OK(SLIP_permute_x(x, n, 1, S));
    int checky = 1;
    if (checky == 1)
    {
        SLIP_info my_check = SLIP_check_solution(A,x,b);
        
        if (my_check == SLIP_OK)
        {
            printf("\nSolution is verified to be exact\n");
        }
        else if (my_check == SLIP_INCORRECT)
        {
            printf("\nERROR! Solution is wrong.\n");
        }
        else
        {
            printf("\nHere");
            return 0;
        }
    }
    
    SLIP_scale_x(x, A, b);

    
    //--------------------------------------------------------------------------
    // Turn Solution into Double (if desired)
    //--------------------------------------------------------------------------
    /*double* x_doub = (double*) SLIP_malloc(n* sizeof(double));
    for (int j = 0; j < n; j++)
        x_doub[j] = mpq_get_d(x[j][0]);    
    
    SLIP_free(x_doub);*/
    
    //--------------------------------------------------------------------------
    // Output & Timing Stats
    //--------------------------------------------------------------------------
    
    double t_col = (double) (end_col-start_col)/CLOCKS_PER_SEC;
    double t_sym = (double) (end_sym-start_sym)/CLOCKS_PER_SEC;
    double t_factor = (double) (end_factor - start_factor) / CLOCKS_PER_SEC;
    double t_solve =  (double) (end_solve - start_solve) / CLOCKS_PER_SEC;

    printf("\nNumber of L nonzeros: \t\t\t%d",
        (L->nz) );
    printf("\nSymmetry Check time: \t\t\t%lf", t_sym);
    printf("\nSymbolic analysis time: \t\t%lf", t_col);
    printf("\nLeft-Looking Chol Factorization time: \t%lf", t_factor);
    printf("\nFB Substitution time: \t\t\t%lf\n\n", t_solve);

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE; 
}
    
