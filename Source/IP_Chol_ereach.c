//------------------------------------------------------------------------------
// IP_Chol/IP_Chol_ereach: Compute reach of an elimint64_t*ation tree
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* Purpose: This function computes the reach of the kth row of A onto the graph of L usint64_t*g the 
   elimint64_t*ation tree. It fint64_t*ds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
   
int64_t IP_Chol_ereach 
(
    SLIP_matrix *A,    // Matrix to be analyzed
    int64_t k,          // Node to start at
    int64_t*  parent,    // ELimint64_t*ation Tree
    int64_t * s,         // Containt64_t*s the nonzero pattern int64_t* s[top..n-1]
    int64_t * w          // Workspace array
)
{
    int64_t i, p, n, len, top ;
    if (!A || !parent || !s || !w) return (-1) ;   /* check int64_t*puts */
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
    return (top) ;                  /* s [top..n-1] containt64_t*s pattern of L(k,:)*/
}
