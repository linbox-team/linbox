// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================

#include "gmp++_int.h"


//-------------------------------------------------- operator -
Integer& Integer::subin(Integer& res, const Integer& n) 
{
  if (iszero(n)) return res;
  if (iszero(res)) return res = - n;
  mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
  return res;
}
Integer& Integer::subin(Integer& res, const long n) 
{
  if (iszero(n)) return res;
  if (iszero(res)) return res = - n;
  int sgn = GMP__SGN(n); 
  if (sgn >0) mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, n);
  else mpz_add_ui((mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, -n);
  return res;
}
Integer& Integer::subin(Integer& res, const unsigned long n) 
{
  if (iszero(n)) return res;
  if (iszero(res)) return res = - n;
  mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
  return res;
}

Integer& Integer::sub(Integer& res, const Integer& n1, const Integer& n2)
{
  if (iszero(n1)) return res = - n2;
  if (iszero(n2)) return res = n1;
  mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
  return res;
}
Integer& Integer::sub(Integer& res, const Integer& n1, const long n2)
{
  if (iszero(n1)) return res = - n2;
  if (iszero(n2)) return res = n1;
  int sgn = GMP__SGN(n2); 
  if (sgn >0) mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
  else mpz_add_ui((mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, -n2);
  return res;
}
Integer& Integer::sub(Integer& res, const Integer& n1, const unsigned long n2)
{
  if (iszero(n1)) return res = - n2;
  if (iszero(n2)) return res = n1;
  mpz_sub_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
  return res;
}

Integer& Integer::neg(Integer& res, const Integer& n)
{
  mpz_neg( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep);
  return res;
}

Integer& Integer::negin(Integer& res)
{
  mpz_neg( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep);
  return res;
}


Integer& Integer::operator -= (const Integer& n)
{
  if (iszero(n)) return *this;
  if (iszero(*this)) return logcpy(-n);
//   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );  
  mpz_sub( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return *this;
}

Integer& Integer::operator -= (const unsigned long l)
{
  if (l==0) return *this;
  if (iszero(*this)) return logcpy(Integer(-l));
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
  mpz_sub_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  return *this;
}

Integer& Integer::operator -= (const long l)
{
  if (l==0) return *this;
  if (iszero(*this)) return logcpy(Integer(-l));
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
  int sgn = GMP__SGN(l);
  if (sgn >0) mpz_sub_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  else mpz_add_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, -l);
  return *this;
}


Integer Integer::operator - (const Integer& n) const
{
  if (iszero(n)) return *this;
  if (iszero(*this)) return -n;
//   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );   
  Integer res;   
  mpz_sub( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return res;
}

Integer Integer::operator - (const unsigned long l) const 
{
  if (l==0) return *this;
  if (iszero(*this)) return Integer(-l);
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
  Integer res;   
  mpz_sub_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l);
  return res;
}

Integer Integer::operator - (const long l) const 
{
  if (l==0) return *this;
  if (iszero(*this)) return Integer(-l);
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
  Integer res;   
  int sgn = GMP__SGN(l);
  if (sgn >0) mpz_sub_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l);
  else mpz_add_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, -l);
  return res;
}
