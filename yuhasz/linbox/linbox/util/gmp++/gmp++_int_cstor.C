#ifndef __GMPplusplus_CSTOR_C__
#define __GMPplusplus_CSTOR_C__
// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================
#include <iostream>
#include "gmp++_int.h"


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
  /* I don't understand what cpt is.  Suppose we just do the else clause?
  if (cpt->decr() ==0) { 
      mpz_clear((mpz_ptr)&gmp_rep) ; 
      delete cpt ;
      mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
  } else
  -bds */
      mpz_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
  return *this ;
}

#endif 

