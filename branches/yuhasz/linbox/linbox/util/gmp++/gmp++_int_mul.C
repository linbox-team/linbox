// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================

#include "gmp++_int.h"


//-------------------------------------------------- operator *
Integer& Integer::mulin(Integer& res, const Integer& n) 
{
  if (iszero(n)) return res = Integer::zero;
  if (iszero(res)) return res;
  mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
  return res;
}
Integer& Integer::mulin(Integer& res, const long n) 
{
  if (iszero(n)) return res = Integer::zero;
  if (iszero(res)) return res;
//   int sgn = GMP__SGN(n); 
  int sgn = GMP__SGN(n); 
  mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, GMP__ABS(n));
//   if (sgn <0) res.gmp_rep.size = -res.gmp_rep.size;
  if (sgn <0) return res = -res;
  return res;
}
Integer& Integer::mulin(Integer& res, const unsigned long n) 
{
  if (iszero(n)) return res = Integer::zero;
  if (iszero(res)) return res;
  mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
  return res;
}

Integer& Integer::mul(Integer& res, const Integer& n1, const Integer& n2)
{
  if (iszero(n1)) return res = Integer::zero;
  if (iszero(n2)) return res = Integer::zero;
  mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
  return res;
}
Integer& Integer::mul(Integer& res, const Integer& n1, const long n2)
{
  if (iszero(n1)) return res = Integer::zero;
  if (iszero(n2)) return res = Integer::zero;
  int sgn = GMP__SGN(n2); 
  mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, GMP__ABS(n2));
//   if (sgn <0) res.gmp_rep.size = -res.gmp_rep.size;
  if (sgn <0) return res = -res;
  return res;
}
Integer& Integer::mul(Integer& res, const Integer& n1, const unsigned long n2)
{
  if (iszero(n1)) return res = Integer::zero;
  if (iszero(n2)) return res = Integer::zero;
  mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
  return res;
}

Integer& Integer::axpy(Integer& res, const Integer& a, const Integer& x, const Integer& b)
{
    if (&res == &b) return Integer::axpyin(res,a,x);
    if (iszero(a) || iszero(x)) return res = b;
    mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&a.gmp_rep, (mpz_ptr)&x.gmp_rep);
    mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&b.gmp_rep);
    return res;
}

Integer& Integer::axpyin(Integer& res, const Integer& a, const Integer& x)
{
    if (iszero(a) || iszero(x)) return res;
    Rep gmp_res; mpz_init((mpz_ptr)&gmp_res);
    mpz_mul( (mpz_ptr)&gmp_res, (mpz_ptr)&a.gmp_rep, (mpz_ptr)&x.gmp_rep);
    mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&gmp_res);
    mpz_clear((mpz_ptr)&gmp_res);
    return res;
}


Integer& Integer::axmy(Integer& res, const Integer& a, const Integer& x, const Integer& b)
{
    if (&res == &b) return Integer::axmyin(res,a,x);
    if (iszero(a) || iszero(x)) return Integer::neg(res,b);
    mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&a.gmp_rep, (mpz_ptr)&x.gmp_rep);
    mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&b.gmp_rep);
    return res;
}
Integer& Integer::axmyin(Integer& res, const Integer& a, const Integer& x)
{
    if (iszero(a) || iszero(x)) return res;
    Rep gmp_res; mpz_init((mpz_ptr)&gmp_res);
    mpz_mul( (mpz_ptr)&gmp_res, (mpz_ptr)&a.gmp_rep, (mpz_ptr)&x.gmp_rep);
    mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&gmp_res);
    mpz_clear((mpz_ptr)&gmp_res);
    return res;
}



Integer& Integer::operator *= (const Integer& n)
{
  if (iszero(n)) return *this = Integer::zero;
  if (iszero(*this)) return *this;
//   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );   
  Integer res; 
  mpz_mul( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return *this = res;
}

Integer& Integer::operator *= (const unsigned long l)
{
  if (l==0) return *this = Integer::zero;
  if (iszero(*this)) return *this;
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
  mpz_mul_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  return *this;
}

Integer& Integer::operator *= (const long l)
{
  if (l==0) return *this =Integer::zero;
  if (iszero(*this)) return *this;
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
  int sgn = GMP__SGN(l);
  mpz_mul_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, GMP__ABS(l));
  if (sgn <0) mpz_neg( (mpz_ptr)&gmp_rep, (mpz_ptr)&(gmp_rep) );
  return *this;
}


Integer Integer::operator * (const Integer& n) const
{
  if (iszero(n)) return Integer::zero;
  if (iszero(*this)) return Integer::zero;
//   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );   
  Integer res;
  mpz_mul( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return res;
}

Integer Integer::operator * (const unsigned long l) const 
{
  if (l==0) return Integer::zero;
  if (iszero(*this)) return Integer::zero;
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
  Integer res;
  mpz_mul_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l);
  return res;
}

Integer Integer::operator * (const long l) const 
{
  if (l==0) return Integer::zero;
  if (iszero(*this)) return Integer::zero;
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
  Integer res;   
  int sgn = GMP__SGN(l);
  mpz_mul_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, GMP__ABS(l));
//   if (sgn <0) (res.gmp_rep).size = -(res.gmp_rep).size;
//   return Integer((res.gmp_rep));
  if (sgn <0) mpz_neg( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&(res.gmp_rep) );
  return res;
}
