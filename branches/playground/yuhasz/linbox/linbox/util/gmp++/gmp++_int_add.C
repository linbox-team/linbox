// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================

#include "gmp++_int.h"

//-------------------------------------------------- operator +
Integer& Integer::addin(Integer& res, const Integer& n) 
{
  if (iszero(n)) return res;
  if (iszero(res)) return res = n;
  mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
  return res;
}
Integer& Integer::addin(Integer& res, const long n) 
{
  if (iszero(n)) return res;
  if (iszero(res)) return res = n;
  int sgn = GMP__SGN(n); 
  if (sgn >0) mpz_add_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, n);
  else mpz_sub_ui((mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, -n);
  return res;
}
Integer& Integer::addin(Integer& res, const unsigned long n) 
{
  if (iszero(n)) return res;
  if (iszero(res)) return res = n;
  mpz_add_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
  return res;
}

Integer& Integer::add(Integer& res, const Integer& n1, const Integer& n2)
{
  if (iszero(n1)) return res = n2;
  if (iszero(n2)) return res = n1;
  mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
  return res;
}
Integer& Integer::add(Integer& res, const Integer& n1, const long n2)
{
  if (iszero(n1)) return res = n2;
  if (iszero(n2)) return res = n1;
  int sgn = GMP__SGN(n2); 
  if (sgn >0) mpz_add_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
  else mpz_sub_ui((mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, -n2);
  return res;
}
Integer& Integer::add(Integer& res, const Integer& n1, const unsigned long n2)
{
  if (iszero(n1)) return res = n2;
  if (iszero(n2)) return res = n1;
  mpz_add_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
  return res;
}


Integer& Integer::operator += (const Integer& n)
{
  if (iszero(n)) return *this;
  if (iszero(*this)) return logcpy(n);
  mpz_add( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return *this;
}

Integer& Integer::operator += (const unsigned long l)
{
  if (l==0) return *this;
  if (iszero(*this)) return logcpy(Integer(l));
  mpz_add_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  return *this;
}

Integer& Integer::operator += (const long l)
{
  if (l==0) return *this;
  if (iszero(*this)) return logcpy(Integer(l));
  int sgn = GMP__SGN(l);
  if (sgn >0) mpz_add_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  else mpz_sub_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, -l);
  return *this;
}


Integer Integer::operator + (const Integer& n) const
{
  if (iszero(n)) return *this;
  if (iszero(*this)) return n;
  Integer res;   
  mpz_add( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return res;
}

Integer Integer::operator + (const unsigned long l) const 
{
  if (l==0) return *this;
  if (iszero(*this)) return Integer(l);
  Integer res;   
  mpz_add_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l);
  return res;
}

Integer Integer::operator + (const long l) const 
{
  if (l==0) return *this;
  if (iszero(*this)) return Integer(l);
  Integer res;   
  int sgn = GMP__SGN(l);
  if (sgn >0) mpz_add_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l);
  else mpz_sub_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, -l);
  return res;
}
