/* Include the Integer-preserving Cholesky routines */

# include "IP-Chol.h"
# include <string.h>

/* Purpose: x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse 
   Arguments:
   A: Input matrix
   j: column of A
   beta: scale
   w: workspace vector
   x: x  = x+ beta A(;,j)
   mark: location of A 
   C: Set in C
   nz: number of nonzeros*/
int SLIP_scatter (SLIP_mat *A, int j, mpz_t beta, int *w, mpz_t *x, int mark,
    SLIP_mat *C, int nz)
{
	int i, p;
	mpz_t temp; mpz_init(temp);
	for (p = A->p [j] ; p < A->p [j+1] ; p++)
	{
		i = A->i [p] ;							// A(i,j) is nonzero 
		if (w [i] < mark)
		{
			w [i] = mark ;						// i is new entry in column j 
			C->i [nz++] = i ;					// add i to pattern of C(:,j) 
			mpz_mul(x[i], beta, A->x[p]);
		}
		else									// x[i] = x[i] + beta*A->x[p]
			mpz_addmul(x[i], beta, A->x[p]);
    }
	mpz_clear(temp);
	return (nz) ;
}

/* Purpose: C = alpha*A + beta*B
   Arguments:
   A: left matrix
   B: right matrix
   alpha: scale of A
   beta: scale of B */
SLIP_mat* SLIP_add (SLIP_mat *A, SLIP_mat *B, mpz_t alpha, mpz_t beta)
{
	/* Allocate memory */
	int p, j, nz = 0, m, n;
	m = A->m ; n = B->n; 
	SLIP_mat *C = new SLIP_mat;
	mpz_t* x = SLIP_initialize_mpz_array(m);
	int* w = SLIP_initialize_int_array(m);
	SLIP_mat_alloc (n, m, A->p[n]+B->p[n], C);
	/* Add C = A+B */
	for (j = 0 ; j < n ; j++)
	{
		C->p [j] = nz ;                   						/* column j of C starts here */
		nz = SLIP_scatter (A, j, alpha, w, x, j+1, C, nz) ;   	/* alpha*A(:,j)*/
		nz = SLIP_scatter (B, j, beta, w, x, j+1, C, nz) ;    	/* beta*B(:,j) */
		for (p = C->p [j] ; p < nz ; p++) mpz_set(C->x[p], x[C->i[p]]);
	}
	C->p [n] = nz ; C->nz = nz;          						/* Finalize the last column of C */
	/* Free memory */
	SLIP_delete_mpz_array(x, m);
	delete[] w;
	SLIP_mat_collapse(C);
	return C;
}

/* Purpose: C = A*B 
   Arguments:
   A: left matrix
   B: right matrix */
SLIP_mat* SLIP_multiply (SLIP_mat *A, SLIP_mat *B)
{
	/* Allocate memory */
	int p, j, nz = 0, m, n;
	m = A->m; n = B->n;
	SLIP_mat *C = new SLIP_mat;
	int* w = SLIP_initialize_int_array(m);
	mpz_t* x = SLIP_initialize_mpz_array(m);
	SLIP_mat_alloc (n, m, A->nz + B->nz, C);
	/* Multiply A*B */
	for (j = 0 ; j < n ; j++)
	{
		if (nz + m > C->nzmax )				// Does C need more room?
			SLIP_mat_realloc(C);
		C->p [j] = nz ;                   	/* column j of C starts here */
		for (p = B->p [j] ; p < B->p [j+1] ; p++)
			nz = SLIP_scatter (A, B->i [p], B->x [p], w, x, j+1, C, nz) ;
		for (p = C->p [j] ; p < nz ; p++) 
			mpz_set(C->x[p], x[C->i[p]]);
		C->nz = nz;
	}
	C->p [n] = nz ; C->nz = nz;
	delete[] w;
	SLIP_delete_mpz_array(x,m);
	SLIP_mat_collapse(C);               	/* remove extra space from C */
	return C ;     							/* success; free workspace, return C */
}




int main( int argc, char* argv[] )
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
    
    /* Symmetric ordering of A */
    //option->order = 2;  // No ordering
    option->order = 1;  // AMD


    
    ok = SLIP_LU_Symbolic(A, S, option, b, 1);
    if (ok != SLIP_OK) return 0;
    if (option->print == 1) SLIP_print_options(option);
    
    /* Uncomment this for hacked AMD */
    /*std::chrono::steady_clock::time_point t_s_start
        = std::chrono::steady_clock::now();
        
    mg_mat* B = SLIP_to_mg(A);
    int j = 0;
    
    //double* w = mg_get_weights(B, 9);
    //double* w2 = mg_get_weights2(A->i, A->p, B->x, n, 9);
 
    //int* p = mg_amd(0, B, 5);
    //int* p = amf(1, B);
    //int * p = mg_amf(0, B, 9);
    //int* p = ammf(1, B);
    
    int* p = (int*) SLIP_malloc(n+1, sizeof(int));
    double Control [AMD_CONTROL];               // Declare AMD control
    amd_defaults(Control);                       // Set AMD defaults
    double Info [AMD_INFO];
    Control[ AMD_DENSE] = 10;
    //Control[ MG_AMD_AGGRESSIVE] = 0;
    /* The options parameter determines how we compute w[k].
     *          0) w[k] = sum (A(i,k)*i) 
     *          1) w[k] = sum (A(:,k))
     *          2) w[k] = mean(A(:,k))
     *          3) w[k] = min (A(:,k))
     *          4) w[k] = max (A(:,k))
     *          5) w[k] = A(k,k), if A(k,k) = 0 first nonzero 
     *          6) w[k] = sum (A(i,k)*(k-i)
     *          7) w[k] = sum (A(i,k)) + n*A(k,k)*/
    /*int options = 9;
    int super = 0;          // 1 means using supernode version

    
    int status = MG_AMD_order(n, A->p, A->i, B->x, p, Control, Info, options, super); // Perform AMD
    //std::cout<<"\nHey the status is: "<<status;
        
    //MG_AMD_order(n, A->p, A->i, p, NULL, NULL);
    
    
     for (int k = 0; k < n; k++)
         S->q[k] = p[k];
         
     std::chrono::steady_clock::time_point t_s_end
         = std::chrono::steady_clock::now();
     std::chrono::duration<float> t_mg = std::chrono::duration_cast<std::chrono::duration<float>> (t_s_end - t_s_start);
//     
//     //std::cout<<"\nNumber of dense columns is: "<< Info[ AMD_NDENSE];
     std::cout<<"\nmg time: \t\t\t"<<t_mg.count();*/
    
    //std::cout<<"\nOur ordering is: \n";
    //for (int k = 0; k < n; k++)
    //    std::cout<<p[k] << " ";
    //SLIP_free(p);
    //mg_mat_delete(B);

    //--------------------------------------------------------------------------
    // Set A2 = PAP'
    //--------------------------------------------------------------------------
    
    /*SLIP_mat* A3 = SLIP_transpose_chol(A);
        
    
    for (int i = 0; i < A->nz; i++)
    {
        //A3->i[i] = A->i[i];
        mpz_neg(A3->x[i], A3->x[i]);
    }
    A3->nz = A->nz;
    
    mpz_t temps; mpz_init(temps); mpz_set_ui(temps, 1);
    SLIP_mat* A4 = SLIP_add(A, A3, temps, temps);
    std::cout<<"\nnnz in A: "<<A->nz;
    std::cout<<"\nnnz in A4: "<<A4->nz;
    
    //for (int i = 0; i < A4->nz; i++)
    //    std::cout<< " " << A4->x[i]; */
    
    int k, index;
    int* pinv2 = (int*) SLIP_malloc(n, sizeof(int));
    for (k = 0; k < n; k++)
    {
        index = S->q[k];
        pinv2[index] = k;
    }
    
    SLIP_mat* A2 = Chol_permute_A(A, pinv2, S);
    A2->p[n] = A->nz;
    A2->nz = A->nz;
    
    //--------------------------------------------------------------------------
    // SLIP LU Factorization
    //--------------------------------------------------------------------------
    Sym_chol* S2;
    S2 = (Sym_chol*) SLIP_malloc(1, sizeof(Sym_chol));
    //ok = SLIP_Chol_Factor(A2, L, S2, rhos, pinv, option);
    //ok = SLIP_Chol_Factor2(A2, L, S2, rhos, pinv, option);
    
    ok = SLIP_Chol_Factor3(A2, L, S2, rhos, pinv, option);
    //std::cout<<"\nok is: "<<ok;
    //std::cout<<"\nLnz is: "<<L->nz <<"\n";
    
    L->m = n;
    if (ok != SLIP_OK) return 0;
//     for (k = 0; k <=n; k++)
//         std::cout<<L->p[k] << " " ;
//     std::cout<<"\nL->i is: \n";
//     for (k = 0; k <L->nz; k++)
//         std::cout<<L->i[k] << " " ;
//     std::cout<<"\nL->x is: \n";
//     for (k = 0; k <L->nz; k++)
//         std::cout<<L->x[k] << " " ;
  

//     U = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
//     SLIP_mat_alloc (L->n, L->m, L->nz, U);
//     
//     // Populate U
//     int count = 0;
//     for (int i = 0; i < n; i++)
//     {
//         for (int j = 0; j < n; j++)
//         {
//             for (int k = L->p[j]; k < L->p[j+1]; k++)
//             {
//                 if (L->i[k] == i)
//                 {
//                     U->i[count] = j;
//                     mpz_set(U->x[count], L->x[k]);
//                     count+=1;
//                 }
//             }
//         }
//         U->p[i+1] = count;
//     } 
//     std::cout<<"\ncount is: "<<count;
//     U->nz = L->nz;
//     
//     SLIP_mat* L2 = (SLIP_mat*) SLIP_malloc(1, sizeof(SLIP_mat));
//     SLIP_mat_alloc (L->n, L->m, L->nz, L2);
//     
//     // Populate U
//     count = 0;
//     for (int i = 0; i < n; i++)
//     {
//         for (int j = 0; j < n; j++)
//         {
//             for (int k = U->p[j]; k < U->p[j+1]; k++)
//             {
//                 if (U->i[k] == i)
//                 {
//                     L2->i[count] = j;
//                     mpz_set(L2->x[count], U->x[k]);
//                     count+=1;
//                 }
//             }
//         }
//         L2->p[i+1] = count;
//     } 
//     std::cout<<"\ncount is: "<<count;
//     L2->nz = L->nz;
    
//    
//     
//     SLIP_mat* L2 = SLIP_transpose_chol(U);
//     std::cout<<"\nL->nz and L2 is: "<<L->nz << " " << L2->nz;
//     std::cout<<"\n For 0: \n\n";
//     std::cout<<"\nL2->i[0] = "<<L2->i[0];
//     std::cout<<"\nL->i[0] = "<<L->i[0];
//     std::cout<<"\nL2->x[0] = "<<L2->x[0];
//     std::cout<<"\nL->x[0] = "<<L->x[0];
//     std::cout<<"\n\n";
//     std::cout<<"\nSize of col 1: "<<L->p[1];
//     for (int i = 0; i < L->nz; i++)
//     {
//         if ( L2->i[i] != L->i[i])
//         {
//             std::cout<<"\nHey, i is: "<<i;
//             std::cout<<"\nL2->["<<i<<"] is: "<<L2->i[i];
//             std::cout<<"\nL->["<<i<<"] is: "<<L->i[i];
//         }
//         if ( mpz_cmp(L2->x[i], L->x[i]) != 0)
//         {
//             std::cout<<"\nHi i is: "<<i;
//             std::cout<<"\nL2->x["<<i<<"] is: "<<L2->x[i];
//             std::cout<<"\nL->x["<<i<<"] is: "<<L->x[i];
//         }
//     }
    
    //std::cout<<"\nFirst column is: "<<S->q[0];
    //std::cout<<"\nLu scale is: "<<option->LU_scale;
    //std::cout<<"\nFirst pivot is: "<<rhos[0];
    //std::cout<<"\nn is: "<<n;
    //std::cout<<"\nSize of the pivots are: ";
    //for (int i = 0; i < n; i++)
    //    std::cout<<" " << mpz_sizeinbase(rhos[i],2); 
    
    //--------------------------------------------------------------------------
    // Solve linear system
    //--------------------------------------------------------------------------
    
    //SLIP_mat* A4 = SLIP_multiply(L, U);
    //std::cout<<"\nA4->n, nz is: "<< A4->n << " " << A4->nz;
    //std::cout<<"\nA2->n, nz is: "<< A2->n << " " << A2->nz;
    
    mpq_t** x = SLIP_initialize_mpq_mat(n,1);
      
    // Solve LDU x = b
    //U = SLIP_transpose_chol(L);
    ok = Up_Solve(x, b, rhos, L, pinv2, option, 1);
    //ok = SLIP_Chol_Solve(x, b, rhos, L, U, pinv2, option, 1);
    if (ok != SLIP_OK) return 0;

    //--------------------------------------------------------------------------
    // Soln verification
    //--------------------------------------------------------------------------
    // x = Q x
    SLIP_Permute_x(x, n, 1, S->q);
    if (option->check == 1) 
        check = SLIP_LU_Check(A, x, n, 1, b, S->q, option);
    SLIP_Scale_x(x, n, 1, option);

    
    /*double* x_doub = (double*) SLIP_malloc(n, sizeof(double));
    for (int j = 0; j < n; j++)
        x_doub[j] = mpq_get_d(x[j][0]);
    
    std::cout<<"\nName is: "<<option->mat_name;
    
    std::string out = ".soln";
    
    std::string out_files = option->mat_name + out;
    
    std::cout<<"\nOutput file is: "<<out_files;
    
    std::ofstream output(out_files);
    
    //--------------------------------------------------------------------------
    // Full precision rational arithmetic
    //--------------------------------------------------------------------------

    //std::cout<<"\nx[1] is: "<<x[0][0];
    //for (int i = 0; i < n; i++)
    //{
          
   //     output << x_doub[i] << "\n";
    //}
    
    
    //std::cout<<"\nx is: \n";
    //for (int j = 0; j < n; j++)
//      std::cout<<x_doub[j] << "\n";
    SLIP_free(x_doub);*/
    
    //--------------------------------------------------------------------------
    // Output
    //--------------------------------------------------------------------------
    Chol_Print_Stats(option, check, L->nz, n, 1, x);
    /*double tot_bits = 0;
    for (k = 0; k < L->nz; k++)
    {
        tot_bits += mpz_sizeinbase(L->x[k],2);
    }
    std::cout<<"\nNumber of L nonzeros: "<<L->nz;
    std::cout<<"\nTotal bits allocated:\t "<<tot_bits;
    std::cout<<"\nAvg bit-length:\t\t "<< (double) tot_bits/L->nz << "\n\n";*/

    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    SLIP_free(pinv2);
    SLIP_mat_delete(A2);
    SLIP_delete_mpz_mat(b,n,1);
    SLIP_delete_mpz_array(rhos,n);
    SLIP_delete_mpq_mat(x,n,1);
    SLIP_mat_delete(A);
    SLIP_mat_delete(L);
    //SLIP_mat_delete(U);
    SLIP_delete_options(option); 
    SLIP_free(pool);
    SLIP_delete_col(S); 
    SLIP_free(pinv);
    slip_gmp_finalize ( ) ;  
}
    
