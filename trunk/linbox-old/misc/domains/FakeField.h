/* ExampleField.h                 Linbox                   -bds 3/00       */
/* ---------------------------------------------------------- */
#ifndef LB__FakeField
#define LB__FakeField
#include <stdlib.h>
#include "interfaces.h"
//#include "Field_a.h"

class FakeField: virtual public Field_a<double>
// doubles faking it.
{private:

 public: 
  FakeField(int x = 0) 
  {}
  typedef double element; 

// start with CR1_a interface:
  Integers Z;
  inline char* label()const
  { return "FakeField (wrapped double)"; }
  inline bool areEqual(const element& a, const element& b)const 
  { return a == b;}
  inline bool areEq(const element& a, const element& b)const
  { return &a == &b;}
  inline istream& read(istream& is, element& a)const
  { return is >> a; }
  inline ostream& write(ostream& os, const element& a)const
  { return os << a; }
  inline Integer& cardinality(Integer& r)const
  { r = one(); return Z.negin(r);} // -1 for infinite cardinality.
  inline element& random(element& r, const Integer& n)const
  { return r = (element)std::random();} // beware if n >= p.  fix.
  //inline element&  random(const Integer& n)const
  //{ return random(*new(element), n); }
  inline element&  random(element& r)const
  { return r = std::random();} 
  //inline element& random()const
  //{ return random(*new(element));}
  inline element& random(element& r, const element& key)const
  { return r = std::random();} // key ignored!  Fix.
  inline const element& zero()const
  { const element& x = 0; return x;} // fix
  inline bool isZero(const element& a)const
  { return a==0; } 
  inline element& add(element& r, const element& a, const element& b)const
  { return r = a + b; } 
  inline element& neg(element& r, const element& a)const
  { return r = -a; } 
  inline element& sub(element& r, const element& a, const element& b)const
  { return r = a - b; } 
  inline element& addin(element& r, const element& b)const
  { return r = r + b; } 
  inline element& negin(element& r)const
  { return r = -r; } 
  inline element& subin(element& r, const element& b)const
  { return r = r - b; } 
//inline element add(const element& a, const element& b)const
//{element r; return add(r, a, b); }
//inline element neg(const element& a)const
//{element r; return neg(r, a);}
//inline element sub(const element& a, const element& b)const
//{element r; return sub(r, a, b);}
  inline element& Zprod(element& r, const Integer n, const element& a)const
  { 
    long m = to_long(n);
    return r = m*a;
  }
  inline element& Zprodin(const Integer n, element& a)const
  { return Zprod(a,n,a); }
//inline element Zprod(const Integer n, const element& a)const
//{element r; return Zprod(r, n, a); }
  inline Integer& characteristic(Integer& r)const
  { return r = 0; }
  inline element& axpy
    (element& r, const element& a, const element& x, const element& y)const
  { return r = (a*x + y); }
  inline element& mul(element& r, const element& a, const element& b)const
  { return r = (a*b); }
  inline element& axmy
    (element& r, const element& a, const element& x, const element& y)const
  { return r = (a*x - y); }
  inline element& mulin(element& r, const element& b)const
  { return r = (r*b); }
  inline element& axpyinx(const element& a, element& x, const element& y)const
  { return x = (a*x + y); }
  inline element& axpyiny(const element& a, const element& x, element& y)const
  { return y = (a*x + y); }
  inline element& axmyinx(const element& a, element& x, const element& y)const
  { return x = (a*x - y); }
  inline element& axmyiny(const element& a, const element& x, element& y)const
  { return y = (a*x - y); }
//inline element& mul(const element& a, const element& b)const
//{ element r; return mul(r, a, b);}
  inline const element& one()const
  { const element& x = 1; return x;} // fix
  inline bool isOne(const element& a)const
  { return a==1; }

// basic field ops
  inline element& inv(element& r, const element& a)const
  { return r = 1/a;}
  inline element& div(element& r, const element& a, const element& b)const
  { return r = a/b;}
  inline element& invin(element& r)const
  { return r = 1/r; }
  inline element& divin(element& r, const element& a)const
  { return r = r/a; }
//inline element& div(const element& a, const element& b)const
//{ return mul(a, inv(b)); }
};

#endif
