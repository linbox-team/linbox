// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================
// Description: 

#include "giverror.h"
#include "givinteger.h"
#include "givcompgmp.h"

Integer pow(const Integer& n, const unsigned long p)
{
  if (p == 0) return Integer::one;

//   Integer::Rep (res.gmp_rep)(l*ABS(n.gmp_rep.size));
  Integer Res;
  mpz_pow_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p);
  return Res;
#if 0
  int is_assg = 0;

  Integer res = Integer::one;
  Integer puiss  = n;

  while (p != 0) {
    if (p & 0x1) 
        if (is_assg) 
           res *= puiss;
        else { 
           is_assg = 1; 
           res = puiss; 
        }
    if ((p >>= 1) != 0) puiss = puiss * puiss;
  }
  return res;
#endif
}

Integer pow(const Integer& n, const long l) {
  GIVARO_ASSERT( l>=0, "[Integer::pow], negative exponent");
  if (l < 0)  return Integer::zero;
  return pow(n, (unsigned long) ABS(l) );
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
  GIVARO_ASSERT( e>=0, "[Integer::powmod], negative exponent not implemented");
  if (e < 0)  return Integer::zero;
  return powmod (n, (unsigned long)ABS(e), m);
}


Integer powmod(const Integer& n, const Integer& e, const Integer& m)
{
  if (e == 0) return Integer::one;
  GIVARO_ASSERT( e>=0, "[Integer::powmod], negative exponent not implemented");
  if (e < 0)  return Integer::zero;
  Integer Res;
  mpz_powm( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, (mpz_ptr)&e.gmp_rep, (mpz_ptr)&m.gmp_rep);
  return Res;
}

