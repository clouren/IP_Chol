//------------------------------------------------------------------------------
// IP_Chol/SLIP_Chol: Solve an SPD linear system using left-chol
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

/* Include the Integer-preserving Cholesky routines */

#define FREE_WORKSPACE                  \
{                                       \
    SLIP_delete_sparse(&A);             \
    SLIP_delete_sparse(&L);             \
    SLIP_delete_sparse(&A2);            \
    SLIP_FREE(pinv2);                   \
    SLIP_delete_dense(&b);              \
    SLIP_delete_mpz_array(&rhos,n);     \
    SLIP_delete_mpq_mat(&x,n,1);        \
    SLIP_FREE(option);                  \
    SLIP_FREE(pinv);                    \
    SLIP_delete_LU_analysis(&S);        \
    SLIP_FREE(S2->parent);              \
    SLIP_FREE(S2->cp);                  \
    SLIP_FREE(S2);                      \
    SLIP_finalize();                    \
}                                       \

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
    int n=0, check, ok, j, index, k, nz = 0;
    
    SLIP_LU_analysis* S = NULL;
    mpz_t** b2 = NULL;
    mpz_t* rhos = NULL;
    int* pinv = NULL;
    int* pinv2 = NULL;
    SLIP_sparse* A2 = NULL;
    Sym_chol* S2 = NULL;
    mpq_t** x = NULL;
    SLIP_sparse *A = SLIP_create_sparse();
    SLIP_sparse *L = SLIP_create_sparse();
    SLIP_dense *b = SLIP_create_dense();
    
    SLIP_options *option = SLIP_create_default_options();
    
    char* mat_name = "./ExampleMats/2.mat";// Set demo matrix and RHS name
    char* rhs_name = "./ExampleMats/2.mat.soln";
    int rat = 1;
    // Default options. May be changed in SLIP_LU_config.h
    // Process the command line
    OK(IP_process_command_line(argc, argv, option,
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
    
    int32_t* I_in, *J_in;
    double* x_in;
  
    OK(IP_tripread_double(A, mat_file, &I_in, &J_in, &x_in, &n, &nz, option));
    //------------------------------------------------------------------
    // At this point, we have read in i, j, and x arrays and have 
    // allocated memory for the A matrix. The i & j are stored as 
    // int32_t and x is stored as a double array. We conclude by using the
    // appropriate SLIP_build_* to construct our input matrix A
    //------------------------------------------------------------------
    fclose(mat_file);
    OK(SLIP_build_sparse_trip_double(A, I_in, J_in, x_in, n, nz, option));
    SLIP_FREE(I_in); SLIP_FREE(J_in); SLIP_FREE(x_in);
    A->n = n;
    
    b2 = SLIP_create_mpz_mat(n,1);
    pinv = (int*) SLIP_malloc(n* sizeof(int));
    rhos = SLIP_create_mpz_array(n);
    // Create RHS
    for (int k = 0; k < n; k++)
        OK(SLIP_mpz_set_ui(b2[k][0],1));
    
    OK(SLIP_build_dense_mpz(b, b2, n, 1));
    SLIP_delete_mpz_mat(&b2, n, 1);
    
    //--------------------------------------------------------------------------
    // Perform Column ordering
    //--------------------------------------------------------------------------
    clock_t start_col = clock();
    
    //SLIP_col* S = (SLIP_col*) SLIP_malloc(1,sizeof(SLIP_col));
    S = SLIP_create_LU_analysis(n+1);
    //S->q = (int*) SLIP_malloc(n+1, sizeof(int));
    
    // Symmetric ordering of A 
    //option->order = 2;  // No ordering
     option->order = 1;  // AMD
    //option->order = 0; // COLAMD
        
    //ok = SLIP_LU_Symbolic(A, S, option, b, 1);
    OK(SLIP_LU_analyze(S, A, option));
    //if (ok != SLIP_OK) return 0;
    //if (option->print == 1) SLIP_print_options(option);
    
    clock_t end_col = clock();
    
    //--------------------------------------------------------------------------
    // Permute matrix A, that is set A2 = PAP'
    //--------------------------------------------------------------------------
    pinv2 = (int*) SLIP_malloc(n* sizeof(int));
    for (k = 0; k < n; k++)
    {
        index = S->q[k];
        pinv2[index] = k;
    }
    
    A2 = IP_Chol_permute_A(A, pinv2, S);
    
    //--------------------------------------------------------------------------
    // SLIP Chol Factorization
    //--------------------------------------------------------------------------
    clock_t start_factor = clock();
    
    S2 = (Sym_chol*) SLIP_malloc(1* sizeof(Sym_chol));
    OK(IP_Left_Chol_Factor(A2, L, S2, rhos, pinv, option));
        
    L->m = n;
    clock_t end_factor = clock();
    
    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    clock_t start_solve = clock();
    
    x = SLIP_create_mpq_mat(n, 1);
    OK(IP_Solve(x, b->x, rhos, L, pinv2, option, 1));
    
    clock_t end_solve = clock();
    
    //--------------------------------------------------------------------------
    // Soln verification
    //--------------------------------------------------------------------------
    // x = Q x
    OK(SLIP_permute_x(x, n, 1, S));
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
    
    double t_sym = (double) (end_col-start_col)/CLOCKS_PER_SEC;
    double t_factor = (double) (end_factor - start_factor) / CLOCKS_PER_SEC;
    double t_solve =  (double) (end_solve - start_solve) / CLOCKS_PER_SEC;

    printf("\nNumber of L nonzeros: \t\t\t%d",
        (L->nz) );
    printf("\nSymbolic analysis time: \t\t%lf", t_sym);
    printf("\nLeft-Looking Chol Factorization time: \t%lf", t_factor);
    printf("\nFB Substitution time: \t\t\t%lf\n\n", t_solve);

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE; 
}
    
