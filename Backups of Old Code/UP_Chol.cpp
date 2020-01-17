# include "SLIP_chol.h"
//# include "../Column_Order/MG/mg_amd_commercial.h"

# include "../SLIP_LU/v0.6/Headers/SLIP_LU_config.h"

# include "UP_chol.h"
//# include "cholmod.h"

# include <string.h>

/* Possible Cholesky Journals

Tim says try SIMAX and TOMS
*/




int main( int argc, char* argv[] )
{
    mp_set_memory_functions(slip_gmp_allocate, slip_gmp_reallocate, slip_gmp_free);
    
    //--------------------------------------------------------------------------
    // Declare memory & Process Command Line
    //--------------------------------------------------------------------------
    int n, check, ok, j, index, k, nz = 0;
    SLIP_mat *A, *L;
    void* pool = SLIP_malloc(1,sizeof(SLIP_LU_Options));        // Must use placement new for this
    SLIP_LU_Options *option = new (pool) SLIP_LU_Options;
    // Default options. May be changed in SLIP_LU_config.h
    SLIP_Set_Options_Defaults(option);
    // Process the command line
    if( SLIP_LU_process_command_line(argc, argv, option) != SLIP_OK)
        return 0;

    option->print = 0;
    option->print2 = 0;
    
    //--------------------------------------------------------------------------
    // Allocate memory    
    //--------------------------------------------------------------------------
    // Allocate A L and U
    A = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
    L = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
    
    // Read in the matrix specified by option->mat_name and store it in A
    SLIP_read_matrix_double(A, option->mat_name, option);
    n = A->n;
    mpz_t** b = SLIP_initialize_mpz_mat(n,1);
    int* pinv = (int*) SLIP_malloc(n, sizeof(int));
    mpz_t* rhos = SLIP_initialize_mpz_array(n);
    // Create RHS
    for (int k = 0; k < n; k++)
        mpz_set_ui(b[k][0],1);
    
    //--------------------------------------------------------------------------
    // Perform Column ordering
    //--------------------------------------------------------------------------
    SLIP_col* S = (SLIP_col*) SLIP_malloc(1,sizeof(SLIP_col));
    S->q = (int*) SLIP_malloc(n+1, sizeof(int));
    
    // Symmetric ordering of A 
    //option->order = 2;  // No ordering
     option->order = 1;  // AMD
    //option->order = 0; // COLAMD
        
    ok = SLIP_LU_Symbolic(A, S, option, b, 1);
    if (ok != SLIP_OK) return 0;
    if (option->print == 1) SLIP_print_options(option);


    //--------------------------------------------------------------------------
    // Permute matrix A, that is set A2 = PAP'
    //--------------------------------------------------------------------------
    int* pinv2 = (int*) SLIP_malloc(n, sizeof(int));
    for (k = 0; k < n; k++)
    {
        index = S->q[k];
        pinv2[index] = k;
    }
    
    SLIP_mat* A2 = Chol_permute_A(A, pinv2, S);
    
    //--------------------------------------------------------------------------
    // SLIP LU Factorization
    //--------------------------------------------------------------------------
    Sym_chol* S2;
    S2 = (Sym_chol*) SLIP_malloc(1, sizeof(Sym_chol));
    //S2->parent = (int*) SLIP_malloc(n, sizeof(int));
    //S2->cp = (int*) SLIP_malloc(n+1, sizeof(int));
     
    ok = Up_Chol_Factor(A2, L, S2, rhos, option);
     
    L->m = n;
    if (ok != SLIP_OK) return 0;
    
   
    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    
    mpq_t** x = SLIP_initialize_mpq_mat(n,1);
    ok = Up_Solve(x, b, rhos, L, pinv2, option, 1);
    if (ok != SLIP_OK) return 0;

    //--------------------------------------------------------------------------
    // Soln verification
    //--------------------------------------------------------------------------
    // x = Q x
    SLIP_Permute_x(x, n, 1, S->q);
    if (option->check == 1) 
        check = SLIP_LU_Check(A, x, n, 1, b, S->q, option);
    SLIP_Scale_x(x, n, 1, option);

    
  // Print the solution out to a file if desired
    double* x_doub = (double*) SLIP_malloc(n, sizeof(double));
    for (int j = 0; j < n; j++)
        x_doub[j] = mpq_get_d(x[j][0]);
     FILE* fp;
     char *end = (char*) ".soln";
     char *name1 = (char*) option->mat_name.c_str();
     char* name2 = (char*) calloc( strlen( name1) + strlen(end) +1, sizeof(char));
     strcat(name2, name1);
     strcat(name2, end);
     fp = fopen(name2, "w");
     for (int i = 0; i < A->n; i++)
         fprintf(fp, "%.16f\n", x_doub[i]);
     fclose(fp);

    
    //--------------------------------------------------------------------------
    // Output
    //--------------------------------------------------------------------------
    Chol_Print_Stats(option, check, L->nz, n, 1, x);
    
    
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
     Sym_chol_free(S2);
     SLIP_free(pinv2);
     SLIP_mat_delete(A2);
     SLIP_delete_mpz_mat(b,n,1);
     SLIP_delete_mpz_array(rhos,n);
     SLIP_delete_mpq_mat(x,n,1);
     SLIP_mat_delete(A);
     SLIP_mat_delete(L);
     SLIP_delete_options(option); 
     SLIP_free(pool);
     SLIP_delete_col(S); 
     SLIP_free(pinv);
     slip_gmp_finalize ( ) ; 
}
    
