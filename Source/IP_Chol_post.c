//------------------------------------------------------------------------------
// REF_Chol/IP_Chol_post: Postorder a forest
//------------------------------------------------------------------------------

// REF Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// Texas A&M University.  All Rights Reserved.  See REF_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/REF-Chol.h"

/* Purpose: post order a forest */
int64_t* IP_Chol_post 
(
    int64_t* parent,    // Parent[j] is parent of node j in forest
    int64_t n           // Number of nodes in the forest
)
{
    int64_t j, k = 0, *post, *w, *head, *next, *stack ;
    if (!parent) return (NULL) ;                                // check inputs 
    post = (int64_t*) SLIP_malloc(n* sizeof(int64_t));          // allocate result 
    w = (int64_t*) SLIP_malloc (3*n* sizeof (int64_t)) ;        // get workspace 
    if (!w || !post) return (NULL) ;
    head = w ; next = w + n ; stack = w + 2*n ;
    for (j = 0 ; j < n ; j++) head [j] = -1 ;           // empty linked lists 
    for (j = n-1 ; j >= 0 ; j--)            // traverse nodes in reverse order
    {
        if (parent [j] == -1) continue ;    // j is a root 
        next [j] = head [parent [j]] ;      // add j to list of its parent 
        head [parent [j]] = j ;
    }
    for (j = 0 ; j < n ; j++)
    {
        if (parent [j] != -1) continue ;    // skip j if it is not a root 
        k = IP_Chol_tdfs (j, k, head, next, post, stack) ;
    }
    SLIP_free(w);
    return (post) ;  // success; free w, return post 
}
