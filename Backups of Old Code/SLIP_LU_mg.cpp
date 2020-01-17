# include "../SLIP_LU/v0.6/Headers/SLIP_LU_config.h"
//# include "../Column_Order/MG/mg.h"
# include "../Column_Order/MG/mg_amd_commercial.h"

/* This program will exactly solve the sparse linear system Ax = b by performing
 * the SLIP LU factorization. Please refer to README.txt for information on how
 * to properly use this code
 */   
int main( int argc, char* argv[])
{
    mp_set_memory_functions(slip_gmp_allocate, slip_gmp_reallocate, slip_gmp_free);
    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int n, check, ok;
    SLIP_mat *A, *L, *U;
    void* pool = SLIP_malloc(1,sizeof(SLIP_LU_Options));        // Must use placement new for this
    SLIP_LU_Options *option = new (pool) SLIP_LU_Options;
    // Default options. May be changed in SLIP_LU_config.h
    SLIP_Set_Options_Defaults(option);
    // Process the command line
    if( SLIP_LU_process_command_line(argc, argv, option) != SLIP_OK)
        return 0;

    
    //--------------------------------------------------------------------------
    // Allocate memory    
    //--------------------------------------------------------------------------
    // Allocate A L and U
    A = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
    L = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
    U = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
    // Read in the matrix specified by option->mat_name and store it in A
    ok = SLIP_read_matrix_double(A, option->mat_name,option);
    if (ok != SLIP_OK) return 0;
    n = A->n;
    mpz_t** b = SLIP_initialize_mpz_mat(n,1);
    int* pinv = (int*) SLIP_malloc(n, sizeof(int));
    mpz_t* rhos = SLIP_initialize_mpz_array(n);
    // Create RHS
    //for (int k = 0; k < n; k++)
    //    mpz_set_ui(b[k][0],1);
    SLIP_read_rhs(b, n, option->rhs_name);

    //--------------------------------------------------------------------------
    // Perform Column ordering
    //--------------------------------------------------------------------------
    SLIP_col* S = (SLIP_col*) SLIP_malloc(1,sizeof(SLIP_col));
    S->q = (int*) SLIP_malloc(n+1, sizeof(int));
    // Column ordering using either AMD, COLAMD, UMFPACK or nothing
    ok = SLIP_LU_Symbolic(A, S, option, b, 1);
    if (ok != SLIP_OK) return 0;
    
    /* Perform mg ordering */
    mg_mat* B  = SLIP_to_mg(A);
    mg_mat* BT = mg_transpose(B);
    
    mg_mat* C = mg_multiply(BT, B);
    
    //int* p = mg_amd(3, B, 5);
    //int* p = amf(3, B);
    //int* p = ammf(3, B);
    
    double alpha = 0.8;
    
    //int * p = mg_amf(3, B, 5, alpha);    
    
    
    int* p = (int*) SLIP_malloc(n+1, sizeof(int));
    double Control [AMD_CONTROL];               // Declare AMD control
    amd_defaults(Control);                       // Set AMD defaults
    double Info [AMD_INFO];
    Control[ AMD_DENSE] = 10;
    int super = 0;
    int options = 5;
    
    int status = MG_AMD_order(n, C->p, C->i, C->x, p, Control, Info, options, super); // Perform AMD

    
    for (int k = 0; k < n; k++)
        S->q[k] = p[k];
    SLIP_free(p);
    mg_mat_delete(B);
    
    if (option->print == 1) SLIP_print_options(option);

    //--------------------------------------------------------------------------
    // SLIP LU Factorization
    //--------------------------------------------------------------------------
    ok = SLIP_LU_Factor(A, L, U, S, rhos, pinv, option);
    if (ok != SLIP_OK) return 0;

    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    mpq_t** x = SLIP_initialize_mpq_mat(n,1);
      
    // Solve LDU x = b
    ok = SLIP_Solve(x, b, rhos, L, U, pinv, option, 1);
    if (ok != SLIP_OK) return 0;

    //--------------------------------------------------------------------------
    // Soln verification
    //--------------------------------------------------------------------------
    // x = Q x
    SLIP_Permute_x(x, n, 1, S->q);
    if (option->check == 1) 
        check = SLIP_LU_Check(A, x, n, 1, b, S->q, option);
    SLIP_Scale_x(x, n, 1, option);

    //--------------------------------------------------------------------------
    // Output
    //--------------------------------------------------------------------------
    SLIP_LU_Print_Stats(option, check, L->nz+U->nz-n, n, 1, x);

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    SLIP_delete_mpz_mat(b,n,1);
    SLIP_delete_mpz_array(rhos,n);
    SLIP_delete_mpq_mat(x,n,1);
    SLIP_mat_delete(A);
    SLIP_mat_delete(L);
    SLIP_mat_delete(U);
    SLIP_delete_options(option); 
    SLIP_free(pool);
    SLIP_delete_col(S); 
    SLIP_free(pinv);
    slip_gmp_finalize ( ) ;
}