#ifndef _INTEGER_H_
#define _INTEGER_H_
// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================
// Description: 
// Integer class definition based on Gmp (>V2.0 or 1.3.2)

extern "C" {
#include "gmp.h"
}
#include "giverror.h"
#include <string>
#include <vector>
class IntPrimeDom;

  //------------------------------------------------------ Class Integer
class Integer {

public:
  typedef vector<mp_limb_t> vect_t;
  Integer( const vector<mp_limb_t>& vect_t );
  //--------------------------------------cstors & dstors
/*Neutral is causing a problem
  Integer(const Neutral n);
*/
  Integer(int n = 0);
  Integer(long n);
  Integer(unsigned int n);
  Integer(unsigned long n);
  Integer(double d);
  Integer(const char *s);
  Integer(const Integer& n);
  ~Integer();

  //------------------------------------- predefined null and one
  static const Integer zero;
  static const Integer one;

  // -- Assignment and copy operators
  Integer& operator = (const Integer& n);
  Integer& logcpy(const Integer& n);
  Integer& copy(const Integer& n);
  
  //------------------Equalities and inequalities between integers and longs
  int operator != (const int l) const;
  int operator != (const long l) const;
  
  friend int compare(const Integer& a, const Integer& b);
  friend int absCompare(const Integer& a, const Integer& b);

  int operator > (const int l) const;
  int operator > (const long l) const;
  int operator < (const int l) const;
  int operator < (const long l) const;

  //----------------Elementary arithmetic between Integers & longs
  Integer& operator += (const Integer& n);  
  Integer& operator += (const unsigned long l);  
  Integer& operator += (const long l);  
  Integer  operator + (const Integer& n) const;  
  Integer  operator + (const unsigned long l) const;
  Integer  operator + (const long l) const;

  Integer& operator -= (const Integer& n);  
  Integer& operator -= (const unsigned long l);  
  Integer& operator -= (const long l);  
  Integer  operator - (const Integer& n) const;
  Integer  operator - (const unsigned long l) const;
  Integer  operator - (const long l) const;
  Integer  operator -() const;

  Integer& operator *= (const Integer& n);  
  Integer& operator *= (const unsigned long l);  
  Integer& operator *= (const long l);  
  Integer  operator * (const Integer& n) const;
  Integer  operator * (const unsigned long l) const;
  Integer  operator * (const long l) const;

  // -- Euclidian division of a/b: returns q or r such that
  // - a=b*q + r, with |r| < |b|, a*r >=0
  Integer& operator /= (const Integer& n);  
  Integer& operator /= (const unsigned long l);
  Integer& operator /= (const long l);
  Integer  operator /  (const Integer& n) const;
  Integer  operator /  (const unsigned long l) const;
  Integer  operator /  (const long l) const;

  Integer& operator %= (const Integer& n);  
  Integer& operator %= (const unsigned long l);
  Integer& operator %= (const long l);
  Integer  operator % (const Integer& n) const;
  long  operator % (const unsigned long l) const;
  long  operator % (const long l) const;

  // - Methods
static Integer& addin (Integer& res, const Integer& n);  
static Integer& addin (Integer& res, const long n);  
static Integer& addin (Integer& res, const unsigned long n);  
static Integer& add   (Integer& res, const Integer& n1, const Integer& n2);  
static Integer& add   (Integer& res, const Integer& n1, const long n2);  
static Integer& add   (Integer& res, const Integer& n1, const unsigned long n2);  

static Integer& subin (Integer& res, const Integer& n);  
static Integer& subin (Integer& res, const long n);  
static Integer& subin (Integer& res, const unsigned long n);  
static Integer& sub   (Integer& res, const Integer& n1, const Integer& n2);  
static Integer& sub   (Integer& res, const Integer& n1, const long n2);  
static Integer& sub   (Integer& res, const Integer& n1, const unsigned long n2);  
static Integer& negin (Integer& res);  
static Integer& neg   (Integer& res, const Integer& n);  

static Integer& mulin (Integer& res, const Integer& n);  
static Integer& mulin (Integer& res, const long n);  
static Integer& mulin (Integer& res, const unsigned long n);  
static Integer& mul   (Integer& res, const Integer& n1, const Integer& n2);  
static Integer& mul   (Integer& res, const Integer& n1, const long n2);  
static Integer& mul   (Integer& res, const Integer& n1, const unsigned long n2);  
static Integer& axpy   (Integer& res, const Integer& a, const Integer& x, const Integer& y );  
static Integer& axpyin   (Integer& res, const Integer& a, const Integer& x);  
static Integer& axmy   (Integer& res, const Integer& a, const Integer& x, const Integer& y );  
static Integer& axmyin   (Integer& res, const Integer& a, const Integer& x);  

static Integer& divin (Integer& q, const Integer& n);  
static Integer& divin (Integer& q, const long n);  
static Integer& divin (Integer& q, const unsigned long n);  
static Integer& div   (Integer& q, const Integer& n1, const Integer& n2);  
static Integer& div   (Integer& q, const Integer& n1, const long n2);  
static Integer& div   (Integer& q, const Integer& n1, const unsigned long n2);  
static Integer& divexact  (Integer& q, const Integer& n1, const Integer& n2);  
static Integer  divexact  (const Integer& n1, const Integer& n2);  

static Integer& modin (Integer& r, const Integer& n);  
static Integer& modin (Integer& r, const long n);  
static Integer& modin (Integer& r, const unsigned long n);  
static Integer& mod   (Integer& r, const Integer& n1, const Integer& n2);  
static Integer& mod   (Integer& r, const Integer& n1, const long n2);  
static Integer& mod   (Integer& r, const Integer& n1, const unsigned long n2);  

  // -- return q, the quotient
static Integer& divmod   (Integer& q, Integer& r, const Integer& n1, const Integer& n2);  
static Integer& divmod   (Integer& q, Integer& r, const Integer& n1, const long n2);  
static Integer& divmod   (Integer& q, Integer& r, const Integer& n1, const unsigned long n2);  
  // -- get it for compatibility
  friend Integer&  divide(const Integer& num, const Integer &den, Integer &q, Integer &r)
  { return Integer::divmod(q,r,num,den); }

  
  //------------------------------------- Arithmetic functions
  friend Integer gcd (const Integer& a, const Integer& b);
  friend Integer gcd (const Integer& a, const Integer& b, 
                            Integer& u, Integer& v);
  friend Integer& gcd (Integer& g, const Integer& a, const Integer& b);
  friend Integer& gcd (Integer& g, const Integer& a, const Integer& b, 
                            Integer& u, Integer& v);

  // - return n^l 
  friend Integer pow(const Integer& n, const long l);
  friend Integer pow(const Integer& n, const unsigned long l);
  friend Integer pow(const Integer& n, const int l) { return pow(n, (long)l ); }
  friend Integer pow(const Integer& n, const unsigned int l) { return pow(n, (unsigned long)l ); }

  // - return n^e % m
  friend Integer powmod(const Integer& n, const unsigned long e, const Integer& m);
  friend Integer powmod(const Integer& n, const long e, const Integer& m);
  friend Integer powmod(const Integer& n, const unsigned int e, const Integer& m) { return powmod(n, (unsigned long)e, m); }
  friend Integer powmod(const Integer& n, const int e, const Integer& m)  { return powmod(n, (long)e, m); }
  friend Integer powmod(const Integer& n, const Integer& e, const Integer& m);

    friend Integer fact ( unsigned long l);
  
    friend Integer sqrt(const Integer& p);
    friend Integer sqrt(const Integer& p, Integer& r);
    friend long logp(const Integer& a, const Integer& p) ;

  //-----------------------------------------Miscellaneous
//   friend int odd(const Integer& n);
//   friend int even(const Integer& n);
  friend Integer abs(const Integer& n);
  friend int probab_prime(const Integer& p);
  friend int probab_prime(const Integer& p, int r);
  friend int jacobi(const Integer& u, const Integer& v) ;
  friend int legendre(const Integer& u, const Integer& v) ;


  Integer operator << (unsigned int l) const; // lshift
  Integer operator >> (unsigned int l) const; // rshift
  Integer operator << (unsigned long l) const; // lshift
  Integer operator >> (unsigned long l) const; // rshift

  // - return the size in byte
  friend inline unsigned long length (const Integer& a); 
  friend inline int sign   (const Integer& a);
  friend inline int iszero (const Integer& a);
  friend inline int isone  (const Integer& a);

  // - return the size in word.
  size_t size() const;
  // - return the i-th word of the integer. Word 0 is lowest word.
  unsigned long operator[](size_t i) const; 

  // -- Convert a Integer to a basic C++ type, throw exception if conversion failed 
  friend long   Integer2long  ( const Integer& n);
  friend vect_t& Integer2vector  (vect_t& v, const Integer& n);
  friend double Integer2double( const Integer& n);
  friend string& Integer2string(string&, const Integer&, int base = 10);
  operator long() const ;
  operator string() const ;
  operator double() const ;
  operator vect_t() const ;


        // -- return a random number with sz machine word. To be improved.
    template< class RandIter > static 
    Integer  random(RandIter& g, int sz=1 );
    template< class RandIter > static 
    Integer  nonzerorandom(RandIter& g, int sz=1 );
    template< class RandIter > static 
    Integer& random(RandIter& g, Integer& r, const Integer& size );
    template< class RandIter > static 
    Integer& nonzerorandom(RandIter& g, Integer& r, const Integer& size );
    template< class RandIter > static 
    Integer& random(RandIter& g, Integer& r, long size =1 );
    template< class RandIter > static 
    Integer& nonzerorandom(RandIter& g, Integer& r, long size =1 );

  //----------------------------------------------I/O

  friend istream& operator >> (istream &i, Integer& n);
  friend ostream& operator << (ostream &o, const Integer&n);
  friend ostream& absOutput (ostream &o, const Integer&n);

  ostream& print( ostream& o ) const;
  
protected:

    friend class IntPrimeDom;
    typedef MP_INT Rep;
    Rep gmp_rep;
    int priv_sign() const;


  // -- Creates a new Integer from a size sz and a array of unsigned long d 
  Integer(unsigned long* d, long size);

public:

protected:
  //---------------------------------------------- Annex functions
  //static const Integer LehmerGcd(const Integer& a, const Integer& b);
static Integer GMPGcd(const Integer& a, const Integer& b);
static Integer& GMPGcd(Integer& g, const Integer& a, const Integer& b);
};//----------------------------------------------- End of Class Integer

// -- return the prime part of Q in P:
extern Integer pp( const Integer& P, const Integer& Q );



//------------------------------------------------------ Class IntegerDom
class IntegerDom : public Integer {
public:
    typedef Integer Rep;
    typedef Rep element;


  IntegerDom() : one(1), zero(0) {}
  IntegerDom(const IntegerDom& D) : one(1), zero(0) {}

  int operator==( const IntegerDom& BC) const { return 1;}
  int operator!=( const IntegerDom& BC) const { return 0;}

  // -- Constants: 
  const Integer one;
  const Integer zero;

  // -- assignement 
  Rep& init  ( Rep& a ) const {return a;}
  Rep& init  ( Rep& a, const Rep& b) const { return a = b ; }
  Rep& read( Rep& a, const long i) const { return a = Integer(i) ; }
  Rep& read( Rep& a, const unsigned long i) const { return a = Integer(i) ; }
  Rep& read( Rep& a, const int i) const { return a = Integer(i) ; }
  Rep& read( Rep& a, const unsigned int i) const { return a = Integer(i) ; }

  Rep& assign( Rep& a, const Rep& b) const { return a = b ; }

  // -- access
  const Rep& access(const Rep& a) const { return a; }
  long   Integer2long  ( const Rep& a) const { return ::Integer2long(a); }
  double Integer2double( const Rep& a) const { return ::Integer2double(a); }
  string& Integer2string(string& s, const Integer& n, int base = 10) const { return ::Integer2string(s,n,base); };

  // -- arithmetic operators
//   Rep& mul( Rep& r, const Rep& a, const Rep& b ) const { return r = a * b; }
//   Rep& div( Rep& r, const Rep& a, const Rep& b ) const { return r = a / b; }
//   Rep& mod( Rep& r, const Rep& a, const Rep& b ) const { return r = a % b; }
//   Rep& add( Rep& r, const Rep& a, const Rep& b ) const { return r = a + b; }
//   Rep& sub( Rep& r, const Rep& a, const Rep& b ) const { return r = a - b; }
  Rep& mul( Rep& r, const Rep& a, const Rep& b ) const { return Integer::mul(r,a,b); }
  Rep& div( Rep& r, const Rep& a, const Rep& b ) const { return Integer::div(r,a,b); }
  Rep& mod( Rep& r, const Rep& a, const Rep& b ) const { return Integer::mod(r,a,b); }
  Rep& add( Rep& r, const Rep& a, const Rep& b ) const { return Integer::add(r,a,b); }
  Rep& sub( Rep& r, const Rep& a, const Rep& b ) const { return Integer::sub(r,a,b); }
  Rep& divmod( Rep& q, Rep& r, const Rep& a, const Rep& b ) const 
  { return Integer::divmod(q,r,a,b); }
  Rep& divexact( Rep& q, const Rep& a, const Rep& b ) const { return Integer::divexact(q,a,b); }

  Rep& mulin( Rep& r, const Rep& a) const { return r *= a; }
  Rep& divin( Rep& r, const Rep& a) const { return r /= a; }
  Rep& modin( Rep& r, const Rep& a) const { return r %= a; }
  Rep& addin( Rep& r, const Rep& a) const { return r += a; }
  Rep& subin( Rep& r, const Rep& a) const { return r -= a; }

//   Rep& axpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const 
//   { return r = a * b + c; }
//   Rep& axmy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const 
//   { return r = a * b - c; }
  Rep& axpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const 
  { return Integer::axpy(r,a,b,c); }
  Rep& axmy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const 
  { return Integer::axmy(r,a,b,c); }
  Rep& axpyin( Rep& r, const Rep& a, const Rep& b ) const 
  { return r += a * b; }
  Rep& axmyin( Rep& r, const Rep& a, const Rep& b ) const 
  { return r -= a * b; }

  // -- unary methods
//   Rep& neg( Rep& r, const Rep& a ) const { return r = -a; }
//   Rep& negin( Rep& r ) const { return r = -r; }
  Rep& neg( Rep& r, const Rep& a ) const { return Integer::neg(r,a); }
  Rep& negin( Rep& r ) const { return Integer::negin(r); }

  // -- extended gcd  q = gcd(a,b) = u*a+v*b;
  Rep& gcd( Rep& g, Rep& u, Rep& v, const Rep& a, const Rep& b ) const 
  { return ::gcd(g, a, b, u, v); }
  Rep& gcd( Rep& g, const Rep& a, const Rep& b ) const 
  { return ::gcd(g, a, b); }

  // - return n^l 
  Rep& pow(Rep& r, const Rep& n, const long l) const { return r = ::pow(n, l); }
  Rep& pow(Rep& r, const Rep& n, const unsigned long l) const { return r = ::pow(n, l); }
  Rep& pow(Rep& r, const Rep& n, const int l) const { return r = ::pow(n, l); }
  Rep& pow(Rep& r, const Rep& n, const unsigned int l) const { return r = ::pow(n, l); }

  // - return square root of n  
  Rep& sqrt(Rep& s, const Rep& n) const { return s = ::sqrt(n); }
  Rep& sqrt(Rep& s, Rep& r, const Rep& n) const { return s = ::sqrt(n, r); }
  // - base p logarithm of a
  long logp(const Rep& a, const Rep& p) const { return ::logp(a,p); }

  // - return n^e % m
  Rep& powmod(Rep& r, const Rep& n, const long e, const Rep& m) const
  { return r = ::powmod(n, e, m);}
  Rep& powmod(Rep& r, const Rep& n, const Rep& e, const Rep& m) const
  { return r = ::powmod(n, e, m);}

  // - Misc
  unsigned long length (const Rep& a) const { return ::length(a); }
  int sign   (const Rep& a) const { return ::sign(a); }
  int iszero (const Rep& a) const { return ::iszero(a); }
  int isone  (const Rep& a) const { return ::isone(a); }
  int isequal (const Rep& a, const Rep& b) const { return compare(a,b) ==0;}
  int isnequal(const Rep& a, const Rep& b) const { return compare(a,b) !=0;}
    int isgeq(const Rep& a, const Rep& b) const { return compare(a,b) >= 0;}
    int isleq(const Rep& a, const Rep& b) const { return compare(a,b) <= 0;}
    int isgeq(const long b,const Rep& a ) const { return isgeq(Rep(b),a);}
    int isleq(const long b,const Rep& a ) const { return isleq(Rep(b),a);}
    int isgeq(const Rep& a, const long b) const { return isgeq(a,Rep(b));}
    int isleq(const Rep& a, const long b) const { return isleq(a,Rep(b));}
    int isgt(const Rep& a, const Rep& b) const { return compare(a,b) > 0;}
    int islt(const Rep& a, const Rep& b) const { return compare(a,b) < 0;}
    int isgt(const long b,const Rep& a ) const { return isgt(Rep(b),a);}
    int islt(const long b,const Rep& a ) const { return islt(Rep(b),a);}
    int isgt(const Rep& a, const long b) const { return isgt(a,Rep(b));}
    int islt(const Rep& a, const long b) const { return islt(a,Rep(b));}

    
    template< class RandIter > Rep& random(RandIter& g, Rep& r, long s = 1) const { return Integer::random(g,r,s); }
    template< class RandIter > Rep& random(RandIter& g, Rep& r, const Rep& b) const { return Integer::random(g,r,b); }
    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, long s = 1) const { return Integer::nonzerorandom(g,r,s); };
    template< class RandIter > Rep& nonzerorandom (RandIter& g,Rep& r, const Rep& b) const { return Integer::nonzerorandom(g,r,b); };

  // -- IO
  istream& read ( istream& i ) 
  { char ch;
    i >> ws >> ch; 
    if (ch != 'I') 
      GivError::throw_error(GivBadFormat("IntegerDom::read: bad signature domain"));
    return i;
  }
  ostream& write( ostream& o ) const { return o << 'I'; }
  istream& read ( istream& i, Rep& n) const { return i >> n; }
  ostream& write( ostream& o, const Rep& n) const { return o << n; }
};


#include "givinteger.inl"

#endif //__INT_H_

