// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================

#include "gmp++_int.h"

// returns 1 if a > b, 0 if a == b and -1 otherwise.  
int compare(const Integer &a, const Integer& b) 
{
  if (&a == &b) return 0 ;
  int cmp = mpz_cmp ( (mpz_ptr)&a.gmp_rep, (mpz_ptr)&b.gmp_rep );
  return (cmp <0 ? -1 : (cmp >0 ? 1 : 0));
}

int absCompare(const Integer &a, const Integer &b) 
{
  Integer c = ((a<0)?-a:a), d = ((b<0)?-b:b);
  return mpz_cmp( (mpz_ptr)&(c.gmp_rep), (mpz_ptr)&(d.gmp_rep));
}

int Integer::operator != (const int l) const 
{ return mpz_cmp_si ( (mpz_ptr)&gmp_rep, l ) != 0; }

int Integer::operator != (const long l) const 
{ return mpz_cmp_si ( (mpz_ptr)&gmp_rep, l ) != 0; }

int Integer::operator > (const int l) const 
{ return mpz_cmp_si((mpz_ptr)&gmp_rep, l) > 0; }

int Integer::operator > (const long l) const 
{ return mpz_cmp_si((mpz_ptr)&gmp_rep, l) > 0; }

int Integer::operator < (const int l) const 
{ return mpz_cmp_si((mpz_ptr)&gmp_rep, l) < 0; }

int Integer::operator < (const long l) const 
{ return mpz_cmp_si((mpz_ptr)&gmp_rep, l) < 0; }

