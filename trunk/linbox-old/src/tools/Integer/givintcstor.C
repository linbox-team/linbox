// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================
#include <iostream.h>
#include "givinteger.h"
#include "givcompgmp.h"


//------------------------------------- predefined null and one
const Integer Integer::zero(0UL);
const Integer Integer::one(1UL);


// -- Integer(const char *s)
Integer::Integer(const char *s) 
{
  mpz_init_set_str((mpz_ptr)&gmp_rep, s, 10);
}


Integer& Integer::copy(const Integer &n)
{
  if (this == &n) return *this;
  if (cpt->decr() ==0) { 
      mpz_clear((mpz_ptr)&gmp_rep) ; 
      delete cpt ;
      mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
  } else
      mpz_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
  return *this ;
}

    
