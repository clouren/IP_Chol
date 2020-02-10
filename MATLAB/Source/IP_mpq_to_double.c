//------------------------------------------------------------------------------
// IP_Chol/IP_mpq_to_double: Convert a mpq_t to double
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "IP-Chol_mex.h"

/* Purpose: This function converts mpq array to double
 * NOTE: This induces roundoff error via the final division
*/
void IP_mpq_to_double
(
    double* x_doub,       // double array
    const mpq_t* x_mpq,   // mpq array
    const int32_t n       // size of b
)
{
    SLIP_info status;
    for (int32_t i = 0; i < n; i++)    
    {
        IP_MEX_OK(SLIP_mpq_get_d(&(x_doub[i]), x_mpq[i]));
    }
}
