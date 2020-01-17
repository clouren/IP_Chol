//------------------------------------------------------------------------------
// IP_Chol/IP_Sym_chol_free: Free Sym_chol data struct
//------------------------------------------------------------------------------

// IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
// Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/IP-Chol.h"

/* Purpose: Free the Sym_Chol Structure*/
void IP_Sym_chol_free
(
    Sym_chol* S
)
{
    SLIP_FREE(S->parent);
    SLIP_FREE(S->cp);
    SLIP_FREE(S);
}