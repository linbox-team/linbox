#ifndef _GIVARO_ZPZ64STD_H_ 
#define _GIVARO_ZPZ64STD_H_ 
// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id$
// ==========================================================================
// Description:
//   Arithmetic on Z/pZ, with p a prime number less than 2^64
//   Modulo typedef is a signed long number. In case it was modified
//   then bezout algorithm must be changed (coefficient can be negative).


#include "givbasictype.h"

#include "giverror.h"
#include "givzpz.h"
class Std64 {}; // -- standard arithmetic over 32bits representations.
typedef long long int64;
typedef unsigned long long uint64;
#define GIVARO_MAXUINT64 18446744073709551615ULL

// ==========================================================================
// -- This class implement the standard arithmetic with Modulo elements:
// - The representation of an integer a in Zpz is the value a % p
// ==========================================================================

template<>
class ZpzDom<Std64> {
public:
  // ----- Exported Types and constantes
  typedef int64 Residu_t;                    // - type to store residue
  enum { size_rep = sizeof(Residu_t) };      // - size of the storage type
  // ----- Representation of element of the domain ZpzDom
  typedef int64 Rep;
  typedef int64 element;


  // ----- Representation of vector of the element
  typedef Rep* Array;
  typedef const Rep* constArray;

  // ----- Constantes 
  const Rep zero;
  const Rep one;

  // ----- Constructor 
  ZpzDom();
  ZpzDom( Residu_t p );

  int operator==( const ZpzDom<Std64>& BC) const { return _p == BC._p;}
  int operator!=( const ZpzDom<Std64>& BC) const { return _p != BC._p;}

  // ----- Access to the modulus 
  Residu_t residu() const;
  Residu_t size() const { return _p; }
  Rep access( const Rep a ) const { return a; }


  // ----- Access to the modulus 
  Rep& init( Rep& a ) const;
  void init( const size_t, Array a, constArray b ) const;
  Rep& init( Rep& a, const long i) const ;
  Rep& init( Rep& a, const unsigned long i) const ;
  Rep& init( Rep& a, const long long i) const ;
  Rep& init( Rep& a, const unsigned long long i) const ;
  Rep& init( Rep& a, const int i) const { return init(a,(long)i); }
  Rep& init( Rep& a, const unsigned int i) const { return init(a,(unsigned long)i); }

  // ----- Misc methods
  int areEqual( const  Rep, const Rep) const;
  int areNEqual( const Rep, const Rep) const;
  int iszero( const Rep a ) const;
  int isnzero( const Rep a ) const;
  int isone ( const Rep a ) const;
  size_t length ( const Rep a ) const;
    bool isZero( const Rep a ) const { return iszero(a); }
    bool isOne ( const Rep a ) const { return isone(a); }

  // ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
  Rep& mul (Rep& r, const Rep a, const Rep b) const;
  Rep& div (Rep& r, const Rep a, const Rep b) const;
  Rep& add (Rep& r, const Rep a, const Rep b) const;
  Rep& sub (Rep& r, const Rep a, const Rep b) const;
  Rep& neg (Rep& r, const Rep a) const;
  Rep& inv (Rep& r, const Rep a) const;

  Rep& mulin (Rep& r, const Rep a) const;
  Rep& divin (Rep& r, const Rep a) const;
  Rep& addin (Rep& r, const Rep a) const;
  Rep& subin (Rep& r, const Rep a) const;
  Rep& negin (Rep& r) const;
  Rep& invin (Rep& r) const;

  // ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
  void mul (const size_t sz, Array r, constArray a, constArray b) const;
  void mul (const size_t sz, Array r, constArray a, Rep b) const;

  void div (const size_t sz, Array r, constArray a, constArray b) const;
  void div (const size_t sz, Array r, constArray a, Rep b) const;

  void add (const size_t sz, Array r, constArray a, constArray b) const;
  void add (const size_t sz, Array r, constArray a, Rep b) const;

  void sub (const size_t sz, Array r, constArray a, constArray b) const;
  void sub (const size_t sz, Array r, constArray a, Rep b) const;

  void neg (const size_t sz, Array r, constArray a) const;
  void inv (const size_t sz, Array r, constArray a) const;

  // -- axpy: r <- a * x + y mod p
  void axpy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
  void axpy 
   (const size_t sz, Array r, constArray a, constArray x, constArray c) const;
  // -- axpyin: r <- r + a * x mod p
  void axpyin(Rep& r, const Rep a, const Rep b) const;
  void axpyin 
   (const size_t sz, Array r, constArray a, constArray x) const;

  // -- amxy: r <- c - a * b mod p
  Rep& amxy (Rep& r, const Rep a, const Rep b, const Rep c) const;
 
  // -- axmy: r <- a * x - y mod p
  void axmy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
  void axmy 
   (const size_t sz, Array r, constArray a, constArray x, constArray c) const;
  // -- axmyin: r <- r - a * x mod p
  void axmyin(Rep& r, const Rep a, const Rep b) const;
  void axmyin 
   (const size_t sz, Array r, constArray a, constArray x) const;

  // -- Misc: r <- a mod p
  void assign ( const size_t sz, Array r, constArray a ) const;
/* JGD 26.10.99
  void assign ( Rep& r, const Rep a) const;
  void assign ( Rep& r, const long a ) const;
  void assign ( Rep& r, const unsigned long a ) const;
  void assign ( Rep& r, const int a ) const;
  void assign ( Rep& r, const unsigned int a ) const;
*/
  Rep& assign ( Rep& r, const Rep a) const;
  Rep& assign ( Rep& r, const long a ) const;
  Rep& assign ( Rep& r, const unsigned long a ) const;
  Rep& assign ( Rep& r, const int a ) const;
  Rep& assign ( Rep& r, const unsigned int a ) const;
   // ----- random generators
//     Rep& NONZEROGIVRANDOM(Rep&) const ;
//     Rep& GIVRANDOM(Rep&) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r, long s) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r, const Rep& b) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, long s) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, const Rep& b) const ;


  // <- \sum_i a[i], return 1 if a.size() ==0,
  void reduceadd ( Rep& r, const size_t sz, constArray a ) const; 

  // <- \prod_i a[i], return 1 if a.size() ==0,
  void reducemul ( Rep& r, const size_t sz, constArray a ) const; 

  // <- \sum_i a[i] * b[i] 
  void dotprod ( Rep& r, const size_t sz, constArray a, constArray b ) const; 
  void dotprod ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const; 

  // ----- a -> r: uint64 to double
  void i2d ( const size_t sz, double* r, constArray a ) const; 

  // ----- a -> r % p: double to uint64 % p
  void d2i ( const size_t sz, Array r, const double* a ) const; 

  // --- IO methods
  istream& read ( istream& s );
  ostream& write( ostream& s ) const;
  istream& read ( istream& s, Rep& a ) const;
  ostream& write( ostream& s, const Rep a ) const;
  Integer& write(Integer&, const Rep a ) const;

protected:
  // -- based for modular inverse, d = a*u + b*v
//   static const int64 gcdext ( int64& u, int64& v, const int64 a, const int64 b );
  int64& gcdext (int64& d, int64& u, int64& v, const int64 a, const int64 b ) const;

protected:
  // -- data representation of the domain:
  Residu_t _p;

  static void Init();
  static void End();
};


#include "givzpz64std.inl"

#endif
