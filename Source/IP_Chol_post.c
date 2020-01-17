//------------------------------------------------------------------------------
// IP_Chol/IP_Chol_post: Postorder a forest
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* Purpose: post order a forest */
int *IP_Chol_post 
(
    int* parent,    // Parent[j] is parent of node j in forest
    int n           // Number of nodes in the forest
)
{
    int j, k = 0, *post, *w, *head, *next, *stack ;
    if (!parent) return (NULL) ;                               /* check inputs */
    post = (int*) SLIP_malloc(n* sizeof(int));                 /* allocate result */
    w = (int*) SLIP_malloc (3*n* sizeof (int)) ;               /* get workspace */
    if (!w || !post) return (NULL) ;
    head = w ; next = w + n ; stack = w + 2*n ;
    for (j = 0 ; j < n ; j++) head [j] = -1 ;           /* empty linked lists */
    for (j = n-1 ; j >= 0 ; j--)            /* traverse nodes in reverse order*/
    {
        if (parent [j] == -1) continue ;    /* j is a root */
        next [j] = head [parent [j]] ;      /* add j to list of its parent */
        head [parent [j]] = j ;
    }
    for (j = 0 ; j < n ; j++)
    {
        if (parent [j] != -1) continue ;    /* skip j if it is not a root */
        k = IP_Chol_tdfs (j, k, head, next, post, stack) ;
    }
    SLIP_free(w);
    return (post) ;  /* success; free w, return post */
}