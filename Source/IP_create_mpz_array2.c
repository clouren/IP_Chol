//------------------------------------------------------------------------------
// IP_Chol/IP_create_mpz_array2: Create an mpz array
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

# include "../Include/IP-Chol.h"

/* Purpose: This function creates an mpz array of size n and allocates
 * memory for numbers of bit size prec. If the relative size of numbers is 
 * known ahead of time, this is more efficient than the
 * SLIP_create_mpz_array
 */
mpz_t* IP_create_mpz_array2
(
    int32_t n,     // size of the array
    int32_t size   // Relative size of numbers
)
{
    SLIP_info ok;
    // Check input
    if (n <= 0 || size <= 0) {return NULL;}
    // Malloc space
    mpz_t* x = (mpz_t*) SLIP_calloc(n, SIZE_MPZ);
    if (!x) {return NULL;}
    for (int32_t i = 0; i < n; i++)
    {
        // Allocate x[i] for bit-length size
        if (SLIP_mpz_init2(x[i],size) != SLIP_OK)
        {
            // Out of memory
            SLIP_MPZ_SET_NULL(x[i]);
            SLIP_delete_mpz_array(&x, n);
            return NULL;
        }
    }
    return x;
}
