//------------------------------------------------------------------------------
// IP_Chol/IP_Chol_ereach: Compute reach of an elimination tree
//------------------------------------------------------------------------------

// IP Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// and Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* Purpose: This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
   
int64_t IP_Chol_ereach 
(
    SLIP_matrix *A,      // Matrix to be analyzed
    int64_t k,           // Node to start at
    int64_t*  parent,    // ELimination Tree
    int64_t * s,         // Contains the nonzero pattern in s[top..n-1]
    int64_t * w          // Workspace array
)
{
    int64_t i, p, n, len, top ;
    // Check inputs
    SLIP_REQUIRE(A, SLIP_CSC, SLIP_MPZ);
    if (!parent || !s || !w) return (-1) ;
    top = n = A->n ; 
    // Mark node k as visited
    SLIP_MARK (w, k) ;
    
    for (p = A->p [k] ; p < A->p [k+1] ; p++)
    {
        // A(i,k) is nonzero
        i = A->i [p] ;
        if (i > k) continue ;  // only use upper triangular part of A 
        for (len = 0 ; !SLIP_MARKED (w,i) ; i = parent [i]) // traverse up etree
        {
            s [len++] = i ;           // L(k,i) is nonzero 
            SLIP_MARK (w, i) ;        // mark i as visited 
        }
        while (len > 0) s [--top] = s [--len] ; // push path onto stack 
    }
    for (p = top ; p < n ; p++) SLIP_MARK (w, s [p]) ;    // unmark all nodes 
    SLIP_MARK (w, k) ;                // unmark node k 
    return (top) ;                    // s [top..n-1] containt64_t*s pattern of L(k,:)
}
