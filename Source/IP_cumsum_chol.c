//------------------------------------------------------------------------------
// IP_Chol/IP_cumsum_chol: Cumulative sum of a vector
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* Purpose: p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c 
   From Tim Davis SuiteSparse */
int64_t IP_cumsum_chol 
(
    int64_t *p, 
    int64_t *c, 
    int64_t n
)
{
    int64_t i, nz = 0 ;
    if (!p || !c) return (-1) ;     /* check inputs */
    for (i = 0 ; i < n ; i++)
    {
        p [i] = nz ;
        nz += c [i] ;
        c [i] = p [i] ;             /* also copy p[0..n-1] back int64_t*o c[0..n-1]*/
    }
    p [n] = nz ;
    return (nz) ;                   /* return sum (c [0..n-1]) */
}
