// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================

#include "gmp++_int.h"


//-------------------------------------------------- operator /
Integer& Integer::divin(Integer& res, const Integer& n) 
{
//  if (iszero(n)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(res)) return res;
  mpz_tdiv_q( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
  return res;
}
Integer& Integer::divin(Integer& res, const long n) 
{
//  if (n ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(res)) return res;
  int sgn = GMP__SGN(n); 
  mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, GMP__ABS(n));
  if (sgn <0) return res = -res;
  return res;
}
Integer& Integer::divin(Integer& res, const unsigned long n) 
{
//  if (n ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(res)) return res;
  mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
  return res;
}

Integer& Integer::div(Integer& res, const Integer& n1, const Integer& n2)
{
  if (iszero(n1)) return res = Integer::zero;
//  if (iszero(n2)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  mpz_tdiv_q( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
  return res;
}
Integer& Integer::div(Integer& res, const Integer& n1, const long n2)
{
  if (iszero(n1)) return res = Integer::zero;
//  if (iszero(n2)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  int sgn = GMP__SGN(n2); 
  mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, GMP__ABS(n2));
  if (sgn <0) return res = -res;
  return res;
}
Integer& Integer::div(Integer& res, const Integer& n1, const unsigned long n2)
{
  if (iszero(n1)) return res = Integer::zero;
//  if (iszero(n2)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  mpz_tdiv_q_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
  return res;
}

Integer& Integer::divexact  (Integer& q, const Integer& n1, const Integer& n2)
{
  if (iszero(n1)) return q = Integer::zero;
//  if (iszero(n2)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  mpz_divexact( (mpz_ptr)&(q.gmp_rep),
                (mpz_ptr)&(n1.gmp_rep), (mpz_ptr)&(n2.gmp_rep)) ;
  return q;
}

Integer  Integer::divexact  (const Integer& n1, const Integer& n2)
{
  if (iszero(n1)) return Integer::zero;
//  if (iszero(n2)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  Integer q;
  mpz_divexact( (mpz_ptr)&(q.gmp_rep),
                (mpz_ptr)&(n1.gmp_rep), (mpz_ptr)&(n2.gmp_rep)) ;
  return q;
}


Integer& Integer::operator /= (const Integer& n)
{
//  if (iszero(n)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(*this)) return *this;
  mpz_tdiv_q( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return *this;
}

Integer& Integer::operator /= (const unsigned long l)
{
//  if (l ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(*this)) return *this;
  mpz_tdiv_q_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
  return *this;
}

Integer& Integer::operator /= (const long l)
{
//  if (l ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(*this)) return *this;
  int sgn = GMP__SGN(l);
  mpz_tdiv_q_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, GMP__ABS(l));
  if (sgn <0) mpz_neg( (mpz_ptr)&gmp_rep, (mpz_ptr)&(gmp_rep));
  return *this;
}


Integer Integer::operator / (const Integer& n) const
{
//  if (iszero(n)) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(*this)) return Integer::zero;
  Integer res;   
  mpz_tdiv_q( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
  return res;
}

Integer Integer::operator / (const unsigned long l) const 
{
//  if (l ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(*this)) return Integer::zero;
  Integer res;   
  mpz_tdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l);
  return res;
}

Integer Integer::operator / (const long l) const 
{
//  if (l ==0) {
//    GivMathDivZero("[Integer::/]: division by zero");
//  }
  if (iszero(*this)) return Integer::zero;
  Integer res;   
  int sgn = GMP__SGN(l);
  mpz_tdiv_q_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, GMP__ABS(l));
  if (sgn <0) return -res;
  return res;
}

// -- Euclidian division
Integer& Integer::divmod(Integer& q, Integer& r, const Integer &a, const Integer &b)
{
//  if (iszero(b)) {
//    GivMathDivZero("[Integer::divide]: division by zero");
//  }

  mpz_tdiv_qr( (mpz_ptr)&(q.gmp_rep), (mpz_ptr)&(r.gmp_rep),
               (mpz_ptr)&(a.gmp_rep), (mpz_ptr)&(b.gmp_rep));
  return q;
}

Integer& Integer::divmod(Integer& q, Integer& r, const Integer& a, const long b)
{
//  if (iszero(b)) {
//    GivMathDivZero("[Integer::divide]: division by zero");
//  }
  int sgn = GMP__SGN(b);
  mpz_tdiv_qr_ui( (mpz_ptr)&(q.gmp_rep), (mpz_ptr)&(r.gmp_rep),
                  (mpz_ptr)&(a.gmp_rep), GMP__ABS(b));
  if (sgn <0) return q = -q;
  return q;
}

Integer& Integer::divmod(Integer& q, Integer& r, const Integer& a, const unsigned long b)
{
//  if (iszero(b)) {
//    GivMathDivZero("[Integer::divide]: division by zero");
//  }
  mpz_tdiv_qr_ui( (mpz_ptr)&(q.gmp_rep), (mpz_ptr)&(r.gmp_rep),
                  (mpz_ptr)&(a.gmp_rep), b);
  return q;
}
