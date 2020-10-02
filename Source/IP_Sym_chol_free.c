//------------------------------------------------------------------------------
// REF_Chol/IP_Sym_chol_free: Free Sym_chol data struct
//------------------------------------------------------------------------------

// REF Chol: (c) 2020, Chris Lourenco, United States Naval Academy, Erick Moreno-Centeno
// Texas A&M University.  All Rights Reserved.  See REF_Chol/License for the license.

//------------------------------------------------------------------------------

#include "../Include/REF-Chol.h"
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
