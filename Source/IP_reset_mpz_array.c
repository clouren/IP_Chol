//------------------------------------------------------------------------------
// IP_Chol/IP_reset_mpz_array: Reset an mpz_t array
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

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
)
{
    SLIP_info ok;
    // Check input
    if (!x || n <= 0 || top < 0 || !xi) {return SLIP_INCORRECT_INPUT;}
    // Access the nonzero pattern located in xi[top..n-1]
    // and set nonzero x[i] = 0
    for (int32_t i = top; i < n; i++)
    {
        OK(SLIP_mpz_set_ui(x[xi[i]], 0));
    }
    return SLIP_OK;
}
