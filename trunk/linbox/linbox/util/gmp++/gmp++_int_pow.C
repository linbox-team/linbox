// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================
// Description: 

#ifdef HAVE_CONFIG_H
#  include "linbox-config.h"
#endif

#include "gmp++/gmp++_int.h"

Integer pow(const Integer& n, const unsigned long p)
{
  if (p == 0) return Integer::one;

  Integer Res;
//  mpz_pow_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p);
__gmpz_pow_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p);
  return Res;
}

Integer pow(const Integer& n, const long l) {
  if (l < 0)  return Integer::zero;
  return pow(n, (unsigned long) GMP__ABS(l) );
}


Integer powmod(const Integer& n, const unsigned long p, const Integer& m)
{
  if (p == 0) return Integer::one;
  Integer Res;
  mpz_powm_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p, (mpz_ptr)&m.gmp_rep);
  return Res;
}
Integer powmod(const Integer& n, const long e, const Integer& m)
{
  if (e < 0)  return Integer::zero;
  return powmod (n, (unsigned long)GMP__ABS(e), m);
}


Integer powmod(const Integer& n, const Integer& e, const Integer& m)
{
  if (e == 0) return Integer::one;
  if (e < 0)  return Integer::zero;
  Integer Res;
  mpz_powm( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, (mpz_ptr)&e.gmp_rep, (mpz_ptr)&m.gmp_rep);
  return Res;
}

