// Integers.h  -- meets CR1 interface.  Fuller interface yet to be defined
// built on NTL ZZ's              Linbox            -bds 8/99          
//-------------------------------------------------------------------
#ifndef LB__Integers
#define LB__Integers
#include "NTL/ZZ.h"

/* ******************************************** 
Expected custom:

namespace user // to have local ZZ disambiguated from NTL's ZZ.
{
  Integers ZZ;
  Integer r, a, b;
  ...
  ZZ.mul(r, a, b);
}
******************************************** */
/* explore this idea later...
class Integer : public ZZ
{private:
  ZZ value;
 public:
  Integer(long v = 0) { value.init(INIT_VAL, v); }
  Integer(ZZ v) : value(v) {};
  const Integer& operator=(const Integer& expr){value = expr.value;}
}
*/  

class Integers 
// meets archetype CR1_a  (also used in defining CR1_a!)
{public:
  typedef ZZ element;
  typedef ZZ Integer;
  Integers(): Zero(INIT_VAL, 0), One(INIT_VAL, 1) {}
  const Integers& operator=(const Integers& b) {}
 private:
  /*static*/ const element Zero/*(INIT_VAL, 0)*/;
  /*static*/ const element One/*(INIT_VAL, 1)*/;

 public:
// start with Ring_a interface:
  char* label(){return "Default Integers (wrap of NTL ZZ)"; }
  bool areEqual(const element& a, const element& b)const;
  bool areEq(const element& a, const element& b)const;
  istream& read(istream& is, element& a)const;
  ostream& write(ostream& os, const element& a)const;
  const Integer& cardinality(Integer& r)const;
  element& random(element& r, const Integer& n)const;
  element& random(const Integer& n)const;
  element& random(element& r)const; 
  //element& random()const; 
  //element& random(element& r, const element& key)const;  // notion of key and notion of integer limit should coincide!?
  const element& zero()const;   
  bool isZero(const element& a)const;
  element& add(element& r, const element& a, const element& b)const; 
  element& neg(element& r, const element& a)const; 
  element& sub(element& r, const element& a, const element& b)const; 
  element& addin(element& r, const element& b)const;  /// r += b;
  element& negin(element& r)const;  
  element& subin(element& r, const element& b)const;  
  element& add(const element& a, const element& b)const;
  element& neg(const element& a)const; /// return -b;
  element& sub(const element& a, const element& b)const;
  element& Zprod(element& r, const Integer n, const element& a)const;
  element& Zprodin(const Integer n, element& a)const;
  element& Zprod(const Integer n, const element& a)const;
  const Integer& characteristic(Integer& r)const;
  element& axpy
  (element& r, const element& a, const element& x, const element& y)const;
  element& mul(element& r, const element& a, const element& b)const; 
  element& axmy
  (element& r, const element& a, const element& x, const element& y)const;
  element& mulin(element& r, const element& b)const; 
  element& axpyinx(const element& a, element& x, const element& y)const; 
  element& axpyiny(const element& a, const element& x, element& y)const; 
  element& axmyinx(const element& a, element& x, const element& y)const; 
  element& axmyiny(const element& a, const element& x, element& y)const; 
  element& mul(const element& a, const element& b)const; 
  const element& one()const;
  bool isOne(const element& a)const;

/*
  element& reducemul(element& r, const Vect& v)const; 
  // r <-- product of entries in v, possibly efficient.
  // Returns one if v.length is 0.

  dotprod(element& r, const Vect& u, const Vect& v)const;
  // r <-- inner product of u and v
  // i.e. r <-- reduceadd( M.mul(u,v) ), where M is the appropriate module.
*/
};

typedef Integers::element Integer;

#endif
