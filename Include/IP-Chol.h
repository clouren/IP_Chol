// This software package exactly solves a sparse symmetric positive definite (SPD)
// system of linear equations using one of two Integer-Preserving Cholesky factorizations. 
// This code accompanies the paper (submitted to SIAM Journal on Matrix Analysis and
// Applications (SIMAX))

//    "Exact Solution of Sparse Symmetric Positive Definite Linear Systems via 
//     Integer Preserving Cholesky Factorization", C. Lourenco, E. Moreno-Centeno, 
//     T. Davis, under submission, SIMAX.

//    If you use this code, you must first download and install the GMP, 
//    MPFR, SLIP LU, AMD, and COLAMD libraries. 
//   
//   GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/
//
//   SLIP LU, AMD, and COLAMD are distributed along with IP-Chol; however,
//   they may be independently obtained at:
//
//   SLIP LU can be found at:
//              https://github.com/clouren/SLIP_LU
//              http://www.suitesparse.com
//
//   AMD and COLAMD can be found at:
//              http://www.suitesparse.com
//
//  
//


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco, Erick Moreno-Centeno, and Timothy Davis

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    IP-Chol is free software; you can redistribute it and/or modify
//     it under the terms of either:
//
//        * the GNU Lesser General Public License as published by the
//          Free Software Foundation; either version 3 of the License,
//          or (at your option) any later version.
//
//     or
//
//        * the GNU General Public License as published by the Free Software
//          Foundation; either version 2 of the License, or (at your option) any
//          later version.
//
//    or both in parallel, as here.
//
//    See license.txt for license info.
//
// This software is copyright by Christopher Lourenco, Erick
// Moreno-Centeno and Timothy A. Davis. All Rights Reserved.
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// IP-Chol is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    This software package solves the SPD linear system Ax = b exactly. The key
//    property of this package is that it can exactly solve ALL SPD input systems. 
//    The input matrix and right hand side vectors are stored as either integers, double
//    precision numbers, multiple precision floating points (through the mpfr
//    library) or as rational numbers (as a collection of numerators and
//    denominators using the GMP mpq_t data structure). Appropriate routines
//    within the code transform the input into an integral matrix in compressed
//    column form.

//    This package computes the factorization PAP' = LDL'. Note that we store the
//    "functional" form of the factorization by only storing the matrix L. The user
//    is given some freedom to select the permutation matrix P. The
//    recommended default settings select P using the AMD ordering.
//    Alternative strategies allowed to select P include the COLAMD 
//    ordering or no column permutation (P=I).  

//    The factor L is computed via integer preserving operations via
//    integer-preserving Gaussian elimination. The key part of this algorithm
//    is a REF Sparse triangular solve function which exploits sparsity and symmetry to
//    reduce the number of operations that must be performed.

//    Once L is computed, a simplified version of the triangular solve
//    is performed which assumes the vector b is dense. The final solution
//    vector x is gauranteed to be exact. This vector can be output in one of
//    three ways: 1) full precision rational arithmetic (as a sequence of
//    numerators and denominators) using the GMP mpq_t data type, 2) double
//    precision while not exact will produce a solution accurate to machine
//    roundoff unless the size of the associated solution exceeds double
//    precision (i.e., the solution is 10^500 or something), 3) variable
//    precision floating point using the GMP mpfr_t data type. The associated
//    precision is user defined.



#ifndef IP_Chol
#define IP_Chol




// IP-Chol inherits many functions, GMP wrappers, etc from SLIP LU. Indeed, IP-Chol is part
// of a suite of exact, integer-preserving sparse matrix software.
// It also inherits the same set of external header files as SLIP LU which includes:
//      <stdlib.h>
//      <stdio.h>
//      <stdbool.h>
//      <stdint.h>
//      <stdlib.h>
//      <string.h>
//      <time.h>
//      <gmp.h>
//      <mpfr.h>

#include "../SLIP_LU-master/SLIP_LU/Include/SLIP_LU.h"
#include <math.h>

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Default Parameters-----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Current version of the code
#define IP_CHOL_VERSION "1.0.0"
#define IP_CHOL_VERSION_MAJOR 1
#define IP_CHOL_VERSION_MINOR 0
#define IP_CHOL_VERSION_SUB   0

// Note that IP-Chol inherits the same error and ordering codes as SLIP LU, thus they are omitted here

// Free a pointer and set it to NULL.
#define SLIP_FREE(p)                        \
{                                           \
    SLIP_free (p) ;                         \
    (p) = NULL ;                            \
}                                           \

#ifndef FREE_WORKSPACE
#define FREE_WORKSPACE  
#endif

#define SLIP_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SLIP_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SLIP_FLIP(i) (-(i)-2)
#define SLIP_UNFLIP(i) (((i) < 0) ? SLIP_FLIP(i) : (i))
#define SLIP_MARKED(Ap,j) (Ap [j] < 0)
#define SLIP_MARK(Ap,j) { Ap [j] = SLIP_FLIP (Ap [j]) ; }


// Size of mpz_t, mpq_t and mpfr_t values
#define SIZE_MPZ  sizeof(mpz_t)
#define SIZE_MPQ  sizeof(mpq_t)
#define SIZE_MPFR sizeof(mpfr_t)


// Field access macros for MPZ/MPQ/MPFR struct
// (similar definition in gmp-impl.h and mpfr-impl.h)

#define MPZ_SIZ(x)   ((x)->_mp_size)
#define MPZ_PTR(x)   ((x)->_mp_d)
#define MPZ_ALLOC(x) ((x)->_mp_alloc)
#define MPQ_NUM(x)   mpq_numref(x)
#define MPQ_DEN(x)   mpq_denref(x)
#define MPFR_MANT(x) ((x)->_mpfr_d)
#define MPFR_EXP(x)  ((x)->_mpfr_exp)
#define MPFR_PREC(x) ((x)->_mpfr_prec)
#define MPFR_SIGN(x) ((x)->_mpfr_sign)
#define MPFR_REAL_PTR(x) (&((x)->_mpfr_d[-1])) /*re-define but same result*/
/* Invalid exponent value (to track bugs...) */
#define MPFR_EXP_INVALID \
 ((mpfr_exp_t) 1 << (GMP_NUMB_BITS*sizeof(mpfr_exp_t)/sizeof(mp_limb_t)-2))

/* Macros to set the pointer in mpz_t/mpq_t/mpfr_t variable to NULL. It is best
 * practice to call these macros immediately after mpz_t/mpq_t/mpfr_t variable
 * is declared, and before the mp*_init function is called. It would help to
 * prevent error when SLIP_MP*_CLEAR is called before the variable is
 * successfully initialized.
 */

#define SLIP_MPZ_SET_NULL(x)                \
    MPZ_PTR(x) = NULL;                      \
    MPZ_SIZ(x) = 0;                         \
    MPZ_ALLOC(x) = 0;                       \

#define SLIP_MPQ_SET_NULL(x)                \
    MPZ_PTR(MPQ_NUM(x)) = NULL;             \
    MPZ_SIZ(MPQ_NUM(x)) = 0;                \
    MPZ_ALLOC(MPQ_NUM(x)) = 0;              \
    MPZ_PTR(MPQ_DEN(x)) = NULL;             \
    MPZ_SIZ(MPQ_DEN(x)) = 0;                \
    MPZ_ALLOC(MPQ_DEN(x)) = 0;              \

#define SLIP_MPFR_SET_NULL(x)               \
    MPFR_MANT(x) = NULL;                    \
    MPFR_PREC(x) = 0;                       \
    MPFR_SIGN(x) = 1;                       \
    MPFR_EXP(x) = MPFR_EXP_INVALID;         \

/* GMP does not give a mechanism to tell a user when an mpz, mpq, or mpfr
 * item has been cleared; thus, if mp*_clear is called on an object that
 * has already been cleared, gmp will crash. It is also not possible to
 * set a mp*_t = NULL. Thus, this mechanism modifies the internal GMP
 * size of entries to avoid crashing in the case that a mp*_t is cleared
 * multiple times.
 */

#define SLIP_MPZ_CLEAR(x)                   \
{                                           \
    if ((x) != NULL && MPZ_PTR(x) != NULL)  \
    {                                       \
        mpz_clear(x);                       \
        SLIP_MPZ_SET_NULL(x);               \
    }                                       \
}                                           \

#define SLIP_MPQ_CLEAR(x)                   \
{                                           \
    SLIP_MPZ_CLEAR(MPQ_NUM(x));             \
    SLIP_MPZ_CLEAR(MPQ_DEN(x));             \
}                                           \

#define SLIP_MPFR_CLEAR(x)                  \
{                                           \
    if ((x) != NULL && MPFR_MANT(x) != NULL)\
    {                                       \
        mpfr_clear(x);                      \
        SLIP_MPFR_SET_NULL(x);              \
    }                                       \
}                                           \

//#ifndef ok
//#define SLIP_info ok
//#endif

#define OK(method)                      \
{                                       \
    ok = method ;                       \
    if (ok != SLIP_OK)                  \
    {                                   \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}                                       \

#define SLIP_CHECK(method)              \
{                                       \
    ok = method ;                       \
    if (ok != SLIP_OK)                  \
    {                                   \
        return 0 ;                      \
    }                                   \
}                                       \

//------------------------------------------------------------------------------
// SLIP_matrix macros
//------------------------------------------------------------------------------

// These macros simplify the access to entries in a SLIP_matrix.
// The type parameter is one of: mpq, mpz, mpfr, int64, or fp64.

// To access the kth entry in a SLIP_matrix using 1D linear addressing,
// in any matrix kind (CSC, triplet, or dense), in any type:
#define SLIP_1D(A,k,type) ((A)->x.type [k])

// To access the (i,j)th entry in a 2D SLIP_matrix, in any type:
#define SLIP_2D(A,i,j,type) SLIP_1D (A, (i)+(j)*((A)->m), type)


//------------------------------------------------------------------------------
// Sym_chol is the data structure for symbolic analysis in the IP-Cholesky 
// factorizations. It includes row permutation, elimination tree, and column
// pointers.
//------------------------------------------------------------------------------

typedef struct Sym_chol
{
    int* pinv;      // Row permutation
    int* parent;    // Elimination tree for Cholesky
    int* cp;        // Column pointers for Cholesky
    int lnz;        // Number of nonzeros in Cholesky L
} Sym_chol;
    
    
/* Purpose: Free the Sym_chol data structure */
void IP_Sym_chol_free
(
    Sym_chol* S
);
   
/* Purpose: Permute the matrix A and return A2 = PAP */
SLIP_info IP_Chol_permute_A
(
    SLIP_matrix **A2_handle // Output permuted matrix
    SLIP_matrix* A,        // Initial input matrix
    int* pinv,             // Row permutation
    SLIP_LU_analysis* S    // Column permutation
);

/* Purpose: Compute the elimination tree of A */

int* IP_Chol_etree 
(
    SLIP_matrix* A // Input matrix (must be SPD)
);

/* Purpose: This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. This is more efficient than the SLIP_reach function 
   It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
   
int IP_Chol_ereach 
(
    SLIP_matrix *A,    // Matrix to be analyzed
    int k,          // Node to start at
    int* parent,    // ELimination Tree
    int* s,         // Contains the nonzero pattern in s[top..n-1]
    int* w          // Workspace array
);

/* Purpose: Depth-first search and postorder of a tree rooted at node j */

int IP_Chol_tdfs 
(
    int j,      // Root node
    int k,      
    int* head,  // Head of list
    int* next,  // Next node in the list
    int* post,  // Post ordered tree
    int* stack  // Stack of nodes
);

/* Purpose: post order a forest */
int *IP_Chol_post 
(
    int* parent,    // Parent[j] is parent of node j in forest
    int n           // Number of nodes in the forest
);


/* Purpose: consider A(i,j), node j in ith row subtree and return lca(jprev,j) 
   Used to determine Column counts of cholesky factor*/
int IP_Chol_leaf 
(
    int i, 
    int j, 
    int* first, 
    int* maxfirst, 
    int* prevleaf,
    int* ancestor, 
    int* jleaf
);

/* Purpose: Something*/
static void init_ata 
(
    SLIP_matrix *AT, 
    int* post, 
    int *w, 
    int **head, 
    int **next
);

/* Purpose: Something*/
int *IP_Chol_counts 
(
    SLIP_matrix *A, 
    int *parent, 
    int *post, 
    int ata // Parameter if we are doing A or A^T A. Setit as 0
);

/* Purpose: This function performs the symmetric sparse REF triangular solve. for uplooking
 * Cholesky factorization. i.e., 
 * (LD) x = A(1:k-1,k). 
 */
SLIP_info IP_Up_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,        // Output the beginning of nonzero pattern
    SLIP_matrix* L,              // partial L matrix
    SLIP_matrix* A,              // input matrix
    int k,                    // iteration of algorithm
    int* xi,                  // nonzero pattern vector
    int* parent,              // Elimination tree
    int* c,                   // Column pointers
    SLIP_matrix* rhos,              // sequence of pivots
    int* h,                   // history vector
    SLIP_matrix* x                  // solution of system ==> kth column of L and U
);
























/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
int IP_Up_Chol_Factor         // performs the Up looking Cholesky factorization
(
    SLIP_sparse* A,           // matrix to be factored
    SLIP_sparse* L,           // lower triangular matrix
    Sym_chol* S,           // Permutation, elimination tree etc
    mpz_t* rhos,           // sequence of pivots
    SLIP_options* option// command options
);

/* Purpose: This solves the system L'x = b for Cholesky factorization */
int IP_Chol_ltsolve 
(
    SLIP_sparse *L,    // The lower triangular matrix
    mpz_t **x,      // Solution vector
    int numRHS      // Number of RHS vectors
);

/* Purpose: This function solves the linear system LD^(-1)L' x = b.*/
int IP_Solve               //solves the linear system LD^(-1)L' x = b
(
    mpq_t** x,              // rational solution to the system
    mpz_t** b,              // right hand side vector
    mpz_t* rhos,            // sequence of pivots
    SLIP_sparse* L,            // lower triangular matrix
    int* pinv,              // row permutation
    SLIP_options* option,// command options
    int numRHS              // number of RHS vectors
);

/* Purpose: p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c 
   From Tim Davis SuiteSparse */
double IP_cumsum_chol 
(
    int *p, 
    int *c, 
    int n
);


/* Purpose: This function sets C = A' */
SLIP_sparse* IP_transpose 
(
    SLIP_sparse *A     // Matrix to be transposed
);


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//--------------------------Alternate Left looking----------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------


/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
int IP_Pre_Left_Factor         // performs the Up looking Cholesky factorization
(
    SLIP_sparse* A,
    SLIP_sparse* L,              // partial L matrix
    int* xi,                  // nonzero pattern vector
    int* parent,              // Elimination tree
    Sym_chol * S,           // stores guess on nnz and column permutation
    int* c                   // Column pointers
);


/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
int IP_Left_Chol_triangular_solve // performs the sparse REF triangular solve
(
    SLIP_sparse* L,              // partial L matrix
    SLIP_sparse* A,              // input matrix
    int k,                    // iteration of algorithm
    int* xi,                  // nonzero pattern vector
    mpz_t* rhos,              // sequence of pivots
    int* pinv,                // inverse row permutation
    int* h,                   // history vector
    mpz_t* x,                  // solution of system ==> kth column of L and U
    int* parent,
    int* c
);

/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
int IP_Left_Chol_Factor         // performs the SLIP LU factorization
(
    SLIP_sparse* A,           // matrix to be factored
    SLIP_sparse* L,           // lower triangular matrix
    Sym_chol * S,           // stores guess on nnz and column permutation
    mpz_t* rhos,           // sequence of pivots
    int* pinv,             // inverse row permutation
    SLIP_options* option// command options
);

/*
 * Purpose: This function allocates a SLIP LU matrix of size n*m with array size
 * nzmax. This version does not allocate individual the values in x. As a
 * result, it is more memory efficient, but also less user friendly.
 *
 * See also IP_sparse_alloc.
 */
SLIP_info IP_sparse_alloc2
(
    SLIP_sparse* A, // sparse matrix data structure to be allocated
    int32_t n,      // number of columns
    int32_t m,      // number of rows (recall m=n assumed)
    int32_t nzmax   // size of allocated i and x arrays
);

/* 
 * Purpose: This function allocates a SLIP LU matrix of size n*m with array
 * size nzmax. This function initializes each entry in A->x therefore they are
 * immediately ready to be operated on. This is less efficient but more user
 * friendly.
 *
 * See also IP_sparse_alloc2.
 */
SLIP_info IP_sparse_alloc
(
    SLIP_sparse* A,// sparse matrix data structure to be allocated
    int32_t n,     // number of columns
    int32_t m,     // number of rows (recall m=n assumed)
    int32_t nzmax  // size of allocated i and x arrays
);

/* 
 * Purpose: This function initializes an int vector of size n and sets the value
 * equal to -1. This function is used for the history and pivot vectors. 
 */
SLIP_info IP_reset_int_array
(
    int32_t *h,    // int32_t vector to be reset
    int32_t n      // size of the int32_t vector
);

/* 
 * Purpose: This function resets an int32_t vector of size n and sets each term
 * equal to -1 with nonzero pattern given. This is more efficient than resetting
 * each term individually
 */
SLIP_info IP_reset_int_array2
(
    int32_t *h,    // int32_t vector to be reset
    int32_t n,     // size of h
    int32_t top,   // beginning of nonzero pattern
    int32_t *xi    // nonzero pattern
);

/* 
 * Purpose: This function resets an mpz array of size n with the nonzero pattern
 * given. This is more efficient than iterating accross all nonzeros in vector x
 */
SLIP_info IP_reset_mpz_array
(
    mpz_t *x,      // mpz array to be reset
    int32_t n,     // size of x
    int32_t top,   // beginning of the nonzero pattern
    int32_t *xi    // nonzero pattern
);

mpz_t* IP_create_mpz_array2
(
    int32_t n,     // size of the array
    int32_t size   // Relative size of numbers
);

SLIP_info IP_get_column //extract k-th column from A, i.e., x=A(:,k)
(
    mpz_t* x,       // A(:,k)
    SLIP_sparse* A, // input matrix
    int32_t k       // column to extract
);


SLIP_info IP_forward_sub
(
    SLIP_sparse *L,   // lower triangular matrix
    mpz_t **x,        // right hand side matrix of size n*numRHS
    mpz_t *rhos,      // sequence of pivots used in factorization
    int32_t numRHS    // number of columns in x
);

SLIP_info IP_array_mul // multiplies vector x by the determinant of matrix
(
    mpz_t** x,      // matrix to be multiplied
    mpz_t det,      // given determinant of matrix
    int32_t n,      // size of x
    int32_t numRHS  // number of RHS vectors
);

SLIP_info IP_array_div // divides the x vector by the determinant
(
    mpq_t** x2,     // solution of x/det
    mpz_t** x,      // input vector
    mpz_t det,      // given determinant of matrix
    int32_t n,      // size of x and x2 
    int32_t numRHS  // number of rhs vectors
);

/* Purpose: This processes the command line for user specified options */
SLIP_info IP_process_command_line //processes the command line
(
    int32_t argc,           // number of command line arguments
    char* argv[],           // set of command line arguments
    SLIP_options* option,   // struct containing the command options
    char** mat_name,        // Name of the matrix to be read in
    char** rhs_name,        // Name of the RHS vector to be read in
    int32_t *rat            // data type of output solution.
                            // 1: mpz, 2: double, 3: mpfr
);

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
    int *n,
    int *nz,
    SLIP_options* option
);

/* Purpose: Determine error in IP Chol, demo only */

void IP_determine_error
(
    SLIP_info ok;
);

/* Purpose: Determine if the input A is indeed symmetric prior to factorization.
 * There are two options as to how to determine the symmetry. 
 * By setting the input exhaustive = 1, both the nonzero pattern and the values
 * of the nonzero entries are checked for symmetry. If A passes both of these tests,
 * then we can be sure it is indeed fully symmetric.
 * 
 * If exhaustive is set to any other value, only the nonzero pattern of A is checked,
 * thus we cannot gauranteee that the matrix is indeed fully symmetric as the values
 * of the entries is not checked.
 * 
 * On success, 0 is returned. If the matrix is not symmetric, 1 is returned.
 * 
 */

int IP_determine_symmetry
(
    SLIP_sparse* A,
    int exhaustive
);

#endif
