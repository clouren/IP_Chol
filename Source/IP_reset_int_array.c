//------------------------------------------------------------------------------
// IP_Chol/IP_reset_int_array: Reset a workspace int array
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* 
 * Purpose: This function initializes an int vector of size n and sets the value
 * equal to -1. This function is used for the history and pivot vectors. 
 */
SLIP_info IP_reset_int_array
(
    int32_t *h,    // int32_t vector to be reset
    int32_t n      // size of the int32_t vector
)    
{
    // Check input
    if (!h || n <= 0) {return SLIP_INCORRECT_INPUT;}
    // Update h
    for (int32_t i = 0; i < n; i++)
    {
        h[i] = -1;
    }
    return SLIP_OK;
}
