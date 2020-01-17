//------------------------------------------------------------------------------
// IP_Chol/IP_Chol_ereach: Compute reach of an elimination tree
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"


/* Purpose: This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
   
int IP_Chol_ereach 
(
    SLIP_sparse *A,    // Matrix to be analyzed
    int k,          // Node to start at
    int* parent,    // ELimination Tree
    int* s,         // Contains the nonzero pattern in s[top..n-1]
    int* w          // Workspace array
)
{
    int i, p, n, len, top ;
    if (!A || !parent || !s || !w) return (-1) ;   /* check inputs */
    top = n = A->n ; 
    SLIP_MARK (w, k) ;                /* mark node k as visited */
    for (p = A->p [k] ; p < A->p [k+1] ; p++)
    {
        i = A->i [p] ;                /* A(i,k) is nonzero */
        if (i > k) continue ;        /* only use upper triangular part of A */
        for (len = 0 ; !SLIP_MARKED (w,i) ; i = parent [i]) /* traverse up etree*/
        {
            s [len++] = i ;         /* L(k,i) is nonzero */
            SLIP_MARK (w, i) ;        /* mark i as visited */
        }
        while (len > 0) s [--top] = s [--len] ; /* push path onto stack */
    }
    for (p = top ; p < n ; p++) SLIP_MARK (w, s [p]) ;    /* unmark all nodes */
    SLIP_MARK (w, k) ;                /* unmark node k */
    return (top) ;                  /* s [top..n-1] contains pattern of L(k,:)*/
}