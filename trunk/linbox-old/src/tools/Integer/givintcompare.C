// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================

#include "givinteger.h"
#include "givcompgmp.h"

// returns 1 if a > b, 0 if a == b and -1 otherwise.  
int compare(const Integer &a, const Integer& b) 
{
  if (&a == &b) return 0 ;

/* Old version of Gmp
  if (a.size > b.size) return 1;
  if (a.size < b.size) return -1;

  if ((a.size == 0) && (b.size == 0)) return 0 ;
  if (a.size == 0)   return (b.size <0) ? 1 : -1 ;
  if (b.size == 0) return (a.size <0) ? -1 : 1 ;
*/

  int cmp = mpz_cmp ( (mpz_ptr)&a.gmp_rep, (mpz_ptr)&b.gmp_rep );
  return (cmp <0 ? -1 : (cmp >0 ? 1 : 0));
/*
  return (a.size <0 ? -cmp : cmp ) ;
*/
}

int absCompare(const Integer &a, const Integer &b) 
{
    Integer c = ((a<0)?-a:a), d = ((b<0)?-b:b);
  return mpz_cmp( (mpz_ptr)&(c.gmp_rep), (mpz_ptr)&(d.gmp_rep));
//   return mpz_cmpabs( (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep));
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
