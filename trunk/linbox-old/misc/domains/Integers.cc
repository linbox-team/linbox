/* Integers.cc              Linbox            -bds 3/00           */
/*-------------------------------------------------------------------*/
#include "Integers.h"
#include "NTL/ZZ.h" 
// we'll make NTL our default integer type.
typedef Integers I;

  typedef ZZ element; 

  bool I::areEqual(const element& a, const element& b)const 
  { return a == b;}
  bool I::areEq(const element& a, const element& b)const
  { return &a == &b;}
  istream& I::read(istream& is, element& a)const
  { return is >> a; }
  ostream& I::write(ostream& os, const element& a)const
  { return os << a; }
  const Integer& I::cardinality(Integer& r)const
  { r = one(); return negin(r); }// -1 stands for infinite
  element& I::random(element& r, const Integer& n)const
  { RandomBnd(r,n); return r; }
  element& I::random(const Integer& n)const
  { element& r = *new(element); RandomBnd(r, n); return r; }
  element& I::random(element& r)const 
  { const Integer n(INIT_VAL, 1000000); 
    RandomBnd(r, n);/*return r = zero(); */} 
  // what should the default distribution be?
  // element& I::random(element& r, const element& key)const //key = n
  //  removed because it duplicates a signature above.
  const element& I::zero()const   
  { return Zero; }
  bool I::isZero(const element& a)const
  { return a == Zero; }
  element& I::add(element& r, const element& a, const element& b)const
  { return r = a + b; }
  element& I::neg(element& r, const element& a)const
  { return r = - a; }
  element& I::sub(element& r, const element& a, const element& b)const 
  { return r = a - b; }
  element& I::addin(element& r, const element& b)const  /// r += b;
  { return r = r + b; }
  element& I::negin(element& r)const
  { return r = - r; }
  element& I::subin(element& r, const element& b)const
  { return r = r - b; }
  element& I::add(const element& a, const element& b)const
  { element& r = *new(element); return r = a + b; }
  element& I::neg(const element& a)const /// return -b;
  { element& r = *new(element); return r = -a; }
  element& I::sub(const element& a, const element& b)const
  { element& r = *new(element); return r = a - b; }
  element& I::Zprod(element& r, const Integer n, const element& a)const
  { return r = n * a; }
  element& I::Zprodin(const Integer n, element& a)const
  { return a = n * a; }
  element& I::Zprod(const Integer n, const element& a)const
  { element& r = *new(element); return r = n*a; }
  const Integer& I::characteristic(Integer& r)const
  { return r = Zero; }
  element& I::axpy
  (element& r, const element& a, const element& x, const element& y)const
  { return r = a*x + y; }
  element& I::mul(element& r, const element& a, const element& b)const
  { return r = a*b; }
  element& I::axmy
  (element& r, const element& a, const element& x, const element& y)const
  { return r = a*x - y; }
  element& I::mulin(element& r, const element& b)const
  { return r = r*b; }
  element& I::axpyinx(const element& a, element& x, const element& y)const 
  { return x = a*x + y; }
  element& I::axpyiny(const element& a, const element& x, element& y)const 
  { return y = a*x + y; }
  element& I::axmyinx(const element& a, element& x, const element& y)const 
  { return x = a*x - y; }
  element& I::axmyiny(const element& a, const element& x, element& y)const
  { return y = a*x - y; }
  element& I::mul(const element& a, const element& b)const 
  { element& r = *new(element); return r = a*b; }
  const element& I::one()const
  { return One;}
  bool I::isOne(const element& a)const
  { return a == One; }

/*
  element& I::reducemul(element& r, const Vect& v)const; 
  // r <-- product of entries in v, possibly efficient.
  // Returns one if v.length is 0.

  I::dotprod(element& r, const Vect& u, const Vect& v)const;
  // r <-- inner product of u and v
  // i.e. r <-- reduceadd( M.mul(u,v) ), where M is the appropriate module.
*/

