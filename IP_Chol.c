//------------------------------------------------------------------------------
// IP_Chol/IP_Chol: Solve an SPD linear system using left or up chol
//------------------------------------------------------------------------------

// IP Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// and Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// IP_Chol/License for the license.

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
    SLIP_LU_analysis_free(&S, NULL);    \
    SLIP_FREE(S2->parent);              \
    SLIP_FREE(S2->cp);                  \
    SLIP_FREE(S2);                      \
    SLIP_finalize();                    \
}                                       \

#define DEMO_OK(method)                 \
{                                       \
    ok = method ;                       \
    if (ok != IP_Chol_OK)               \
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
    // Prior to using IP-Chol, its environment must be initialized. It is 
    // initialized using the SLIP LU enviroment, the SLIP_initialize function.
    //--------------------------------------------------------------------------
    
    SLIP_initialize();
    
    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int64_t n = 0, check, ok, j, index, k, nz = 0;
   
    SLIP_LU_analysis* S = NULL;
    SLIP_matrix *A = NULL;
    SLIP_matrix *L = NULL;
    SLIP_matrix *b = NULL;
    SLIP_matrix *rhos = NULL;
    int64_t* pinv = NULL;
    int64_t* pinv2 = NULL;
    SLIP_matrix* A2 = NULL;
    Sym_chol* S2 = NULL;
    SLIP_matrix* x = NULL;
    SLIP_matrix* L2 = NULL;
    SLIP_matrix* rhos2 = NULL;
    S2 = (Sym_chol*) SLIP_malloc(1* sizeof(Sym_chol));
    
    // Default options. May be changed in SLIP_LU_config.h
    SLIP_options *option = SLIP_create_default_options();
    
    char* mat_name = "./ExampleMats/872.mat";// Set demo matrix and RHS name
    char* rhs_name = "./ExampleMats/872.mat.soln";
    int64_t rat = 1;
    bool left = false; // Default is up-looking factorization
    
    // Process the command line
    DEMO_OK(IP_process_command_line(argc, argv, &left, option,
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
    
    DEMO_OK(IP_tripread_double(&A, mat_file, option));
    fclose(mat_file);
    n = A->n;
    // For this code, we utilize a vector of all ones as the RHS vector    
    SLIP_matrix_allocate(&b, SLIP_DENSE, SLIP_MPZ, n, 1, n, false, true, option);
    pinv = (int64_t*) SLIP_malloc(n* sizeof(int64_t));
    // Create RHS
    for (int64_t k = 0; k < n; k++)
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
    // Determine if A is indeed symmetric. The user has two options in this case.
    // By setting the second input argument to true, IP_Chol performs an 
    // exhaustive search for symmetry. That is, it checks both the nonzero pattern,
    // and ensures that the corresponding values are identical. Thus, the matrix 
    // can fail this check if it is structually symmetric but not actually symmetric.
    //
    // Alternatively, if the second argument is set to false, IP_Chol checks if 
    // the matrix is structually symmetric but does not verify that the values are
    // indeed identical.
    // --------------------------------------------------------------------------
    clock_t start_sym = clock();
    
    // TODO Handle error conditions here.
    
    //DEMO_OK(IP_determine_symmetry(A, false));    // Determine symmetry just with nonzero pattern
    DEMO_OK(IP_determine_symmetry(A, true));    // Determine symmetry with nonzero pattern and values
        
    
    clock_t end_sym = clock();
    //--------------------------------------------------------------------------
    // Permute matrix A, that is set A2 = PAP'
    //--------------------------------------------------------------------------
    pinv2 = (int64_t*) SLIP_malloc(n* sizeof(int64_t));
    for (k = 0; k < n; k++)
    {
        index = S->q[k];
        pinv2[index] = k;
    }
    
    DEMO_OK( IP_Chol_permute_A(&A2, A, pinv2, S));
    
    // Debugging code to check the matrix if desired
    option->print_level = 2;
    option->check = true;
    //SLIP_matrix_check(A2,option);
    
    //--------------------------------------------------------------------------
    // IP Chol Factorization
    //--------------------------------------------------------------------------
    clock_t start_factor = clock();
    
    DEMO_OK( IP_Chol_Factor( A2, &L, S2, &rhos, left, option));
//    L->m = n;

    clock_t end_factor = clock();
     
    
    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    clock_t start_solve = clock();
    option->check = true;
    
    DEMO_OK( IP_Solve( &x, A2, A, b, rhos, L, pinv2, S, option));
    
    
    clock_t end_solve = clock();
    
    //--------------------------------------------------------------------------
    // Output & Timing Stats
    //--------------------------------------------------------------------------
    
    double t_col = (double) (end_col-start_col)/CLOCKS_PER_SEC;
    double t_sym = (double) (end_sym-start_sym)/CLOCKS_PER_SEC;
    double t_factor = (double) (end_factor - start_factor) / CLOCKS_PER_SEC;
    double t_solve =  (double) (end_solve - start_solve) / CLOCKS_PER_SEC;

    // TODO: inherit the SLIP LU printing functions for error handling etc
    if (left == true)
        printf("\nLeft-looking Factorization stats:");
    else
        printf("\nUp-looking Factorization stats:");
    printf("\nNumber of L nonzeros: \t\t\t%ld",
        (L->p[L->n]) );
    printf("\nSymmetry Check time: \t\t\t%lf", t_sym);
    printf("\nSymbolic analysis time: \t\t%lf", t_col);
    printf("\nIP Chol Factorization time: \t\t%lf", t_factor);
    printf("\nFB Substitution time: \t\t\t%lf\n\n", t_solve);

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    FREE_WORKSPACE;
}
    
