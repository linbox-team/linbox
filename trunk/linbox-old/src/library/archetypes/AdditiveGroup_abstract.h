/* AdditiveGroup.h                 Linbox                   -bds 3/00 7/00    */
/* ---------------------------------------------------------- */
#ifndef __AdditiveGroup_abstract__
#define __AdditiveGroup_abstract__
#include <iostream.h>
#include <LinBox/BasicDomain_abstract.h>
#include <LinBox/integer.h>

namespace linbox{

/** 
AdditiveGroup_a interface - definition for the category AbelianGroup
with add() as the basic binary operation.

zero(), isZero(), 
add(), neg(), sub(),
addin(), negin(), subin(),
and the Set functionality.
*/
class AdditiveGroup_abstract : public BasicDomain_abstract
{public: 

// Start with Set_abstract interface and add these:
  /** 
  zero() is the identity element: zero() + a == a, zero()*a == zero().
  */
  virtual const eltbase& zero() const
  =0;   

  /** 
  isZero(a) // a == 0
  */
  virtual bool isZero( const eltbase& a ) const
  =0;

  /** 
  add(r,a,b) // r = a + b, commutative, associative.
  */
  virtual eltbase& add( eltbase& r, const eltbase& a, const eltbase& b ) const
  =0; 

  /** 
  neg(r,a) // r = -a.  i.e. isZero(add(a,neg(a))).
  */
  virtual eltbase& neg( eltbase& r, const eltbase& a ) const
  =0; 

  /** 
  sub(r,a,b) // r = a + -b.
  */
  virtual eltbase& sub( eltbase& r, const eltbase& a, const eltbase& b )const
  =0; 

// in place forms
  /**
  addin(r,b) // r += b.
  */
  virtual eltbase& addin( eltbase& r, const eltbase& b ) const
  =0; 

  /**
  negin(r) // r = -r.
  */
  virtual eltbase& negin( eltbase& r ) const
  =0;  

  /**
  subin(r,b) // r -= b.
  */
  virtual eltbase& subin( eltbase& r, const eltbase& b ) const
  =0;  

// scalar mul by Integer, (abelian gp is Z_Module)
  /**
  Zmul(r, n, a) // r <- n*a, scalar product
  */
  virtual eltbase& Zmul( eltbase& r, const Integer n, const eltbase& a ) const
  =0;

  /**
  Zmulin(n, a) // a <- n*a, scalar product
  */
  virtual eltbase& Zmulin( const Integer n, eltbase& a ) const
  =0;

  /**
  n = characteristic() if n is least > 0 such that Zmul(n,a) = 0, for all a.
  0 = characteristic() if no such positive n.
  */
  virtual Integer& characteristic( Integer& r ) const
  =0;

}; // AdditiveGroup_abstract

} //namespace linbox
#endif
