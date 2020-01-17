//------------------------------------------------------------------------------
// IP_Chol/IP_reset_int_array2: Reset an int array to 0s
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"


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
)
{
    if (!h || !xi || n <= 0 || top < 0) {return SLIP_INCORRECT_INPUT;}
    // Access the nonzero pattern located in xin[top..n-1]
    // and set Each "nonzero" h[i] = -1 
    for (int32_t i = top; i < n; i++)
    {
        h[xi[i]] = -1;
    }
    return SLIP_OK;
}
