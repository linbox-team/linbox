/* ExampleField.h                 Linbox                   -bds 3/00       */
/* ---------------------------------------------------------- */
#ifndef LB__ExampleField
#define LB__ExampleField
#include <stdlib.h>
#include "interfaces.h"
//#include "Field_a.h"

class ExampleField: virtual public Field_a<int>
// small prime fields with tight representation
{private:
  const int p;
 public: 
  ExampleField (int prime = 3): p(prime){}
  typedef int element; // has to be explicit here?

// start with CR1_a interface:
  inline char* label()const 
  { return "ExampleField (brute force impl of small prime fields)"; } 
  inline bool areEqual(const element& a, const element& b)const 
  {return a == b;}
  inline bool areEq(const element& a, const element& b)const
  {return &a == &b;}
  inline istream& read(istream& is, element& a)const
  {return is >> a; }
  inline ostream& write(ostream& os, const element& a)const
  {return os << a; }
  inline Integer& cardinality(Integer& r)const
  { return r = p;}
  inline element& random(element& r, const Integer& n)const
  { return r = random()%p;} // beware if n >= p.  fix.
  inline element&  random(const Integer& n)const
  { return random(*new(element), n); }
  inline element&  random(element& r)const
  { return r = std::random()%p;} 
  inline element& random()const
  { return random(*new(element));}
  inline element& random(element& r, const element& key)const
  { return r = std::random()%p;} // key ignored!  Fix.
  inline const element& zero()const
  { const element& x = 0; return x;} // fix
  inline bool isZero(const element& a)const
  {return a==0; } 
  inline element& add(element& r, const element& a, const element& b)const
  {r = a + b; if (r >= p) r -= p; return r; } 
  inline element& neg(element& r, const element& a)const
  {return r = -a + p; } 
  inline element& sub(element& r, const element& a, const element& b)const
  {r = a - b; if (r < 0) r += p; return r; } 
  inline element& addin(element& r, const element& b)const
  {r = r + b; if (r >= p) r -= p; return r; } 
  inline element& negin(element& r)const
  {return r = -r + p; } 
  inline element& subin(element& r, const element& b)const
  {r = r - b; if (r < 0) r += p; return r; } 
//inline element add(const element& a, const element& b)const
//{element r; return add(r, a, b); }
//inline element neg(const element& a)const
//{element r; return neg(r, a);}
//inline element sub(const element& a, const element& b)const
//{element r; return sub(r, a, b);}
  inline element& Zprod(element& r, const Integer n, const element& a)const
  {return r = (n*a)%p;}
  inline element& Zprodin(const Integer n, element& a)const
  {return a = (n*a)%p;}
//inline element Zprod(const Integer n, const element& a)const
//{element r; return Zprod(r, n, a); }
  inline Integer& characteristic(Integer& r)const
  { return r = p; }
  inline element& axpy
    (element& r, const element& a, const element& x, const element& y)const
  {return r = (a*x + y)%p; }
  inline element& mul(element& r, const element& a, const element& b)const
  {return r = (a*b)%p; }
  inline element& axmy
    (element& r, const element& a, const element& x, const element& y)const
  { r = a*x - y; return normalized(r); }
  inline element& mulin(element& r, const element& b)const
  {return r = (r*b)%p; }
  inline element& axpyinx(const element& a, element& x, const element& y)const
  {return x = (a*x + y)%p; }
  inline element& axpyiny(const element& a, const element& x, element& y)const
  {return y = (a*x + y)%p; }
  inline element& axmyinx(const element& a, element& x, const element& y)const
  { x = a*x - y; return normalized(x); }
  inline element& axmyiny(const element& a, const element& x, element& y)const
  { y = a*x - y; return normalized(y); }
//inline element& mul(const element& a, const element& b)const
//{element r; return mul(r, a, b);}
  inline const element& one()const
  { const element& x = 1; return x;} // fix
  inline bool isOne(const element& a)const
  { return a==1; }

// basic field ops
  inline element& inv(element& r, const element& a)const
  { return r = inv(a);}
  inline element& div(element& r, const element& a, const element& b)const
  { return mul(r, a, inv(b)); }
  inline element& invin(element& r)const
  { return r = inv(r); }
  inline element& divin(element& r, const element& a)const
  { return mulin(r, inv(a)); }
//inline element& div(const element& a, const element& b)const
//{ return mul(a, inv(b)); }
 private:
  inline element& normalized(element& a)const
  { while (a < 0) a += p;
    a = a % p;
    return a;
  }
  inline element inv(const element& a)const
  { element s; hegcd(s, a, p); 
    return normalized(s);
  }
  inline void hegcd(element& s, const element& a, const element& b)const
  /// return s such that gcd(a,b) = sa + tb, for some t.
  { element r1 = b; 
    element r2 = a;
    element s1 = 0;
    element s2 = 1;
    while (r2 != 1) 
    { element q = r1/r2;
      element r3 = r1 - q*r2; element s3 = s1 - q*s2;
      r1 = r2; r2 = r3; 
      s1 = s2; s2 = s3;
    }
    s = s2;
  }


};

#endif
