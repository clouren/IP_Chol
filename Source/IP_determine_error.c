///------------------------------------------------------------------------------
// IP_Chol/IP_determine_error: Output error code for IP Chol
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

/* Purpose: Determine what error occured when using IP Chol. This is used solely
 * for demo/debugging purposes.
 */

#include "../Include/IP-Chol.h"

void IP_determine_error
(
    SLIP_info ok
)
{
    if (ok == SLIP_OUT_OF_MEMORY)
    {
        printf("\nOut of memory\n");
    }
    else if (ok == SLIP_SINGULAR)
    {
        printf("\nInput matrix is singular OR no diagonal pivot. Please ensure input is SPD\n");
    }
    else if (ok == SLIP_INCORRECT_INPUT)
    {
        printf("\nIncorrect input for a IP Chol Function\n");
    }
}