/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// ==========================================================================
// file: lin_dom_tt.inl
// Description:
// Bugs:
// Time-stamp: <11 Apr 03 17:00:33 Jean-Guillaume.Dumas@imag.fr> 
// ==========================================================================

 // ------------------------- Miscellaneous functions

template<class TT> inline bool OperatorWrapper<TT>::iszero(const Rep& a) const
  { return a == zero ; }

template<class TT> inline bool OperatorWrapper<TT>::isone(const Rep& a) const
  { return a == one ; }

template<class TT> inline bool OperatorWrapper<TT>::isnzero(const Rep& a) const
  { return a != zero ; }

template<class TT> inline bool OperatorWrapper<TT>::isnone(const Rep& a) const
  { return a != one ; }

template<class TT> inline size_t OperatorWrapper<TT>::length(const Rep& ) const
  { return sizeof(TT) ;}

template<class TT> inline     bool OperatorWrapper<TT>::isequal( const Rep& a, const Rep& b) const { return a == b ; }
    
template<class TT> inline     bool OperatorWrapper<TT>::isnequal( const Rep& a, const Rep& b) const { return a != b ; }
    
template<class TT> inline     bool OperatorWrapper<TT>::islt( const Rep& a, const Rep& b) const { return a < b ; }
    
template<class TT> inline     bool OperatorWrapper<TT>::isgt( const Rep& a, const Rep& b) const { return a > b ; }
    
  // ----------- Usefull method :
template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::mul (Rep& r, const Rep&a, const Rep& b) const 
  { return r = a * b; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::mulin (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a) const 
  { return r *= a; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::div (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a, const OperatorWrapper<TT>::Rep&b) const 
  { return r = a / b; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::divin (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a) const 
  { return r /= a; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::add (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a, const OperatorWrapper<TT>::Rep&b) const 
  { return r = a + b; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::addin (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a) const 
  { return r += a; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::sub (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a, const OperatorWrapper<TT>::Rep&b) const 
  { return r = a - b; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::subin (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a) const 
  { return r -= a; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::neg (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a) const 
  { return r = -a; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::negin (OperatorWrapper<TT>::Rep& r) const 
  { return r = -r; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::inv (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a) const 
  { return r = one/a; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::invin (OperatorWrapper<TT>::Rep& r) const 
  { return r = one/r; }

/*
#include "lin_sqrt.h"

template<class TT> inline OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::sqrt (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a) const 
  { return r = ::sqrt(a); }

template<class TT> inline OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::sqrtin (OperatorWrapper<TT>::Rep& r) const 
  { return r = ::sqrt(r); }
*/

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::axpy (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a, const OperatorWrapper<TT>::Rep&b, const OperatorWrapper<TT>::Rep&c) const 
  { return r = a * b + c; }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::axpyin (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a, const OperatorWrapper<TT>::Rep&b) const { return r = r + a*b ; }

  // -- axmyin: r <- r - a * b mod p
template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::axmyin (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a, const OperatorWrapper<TT>::Rep&b) const { return r -= a*b;  }

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::axmy (OperatorWrapper<TT>::Rep& r, const OperatorWrapper<TT>::Rep&a, const OperatorWrapper<TT>::Rep&b, const OperatorWrapper<TT>::Rep&c)  const { return r = a * b - c; }


template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::assign( Rep& r, const Rep&residu ) const {
    return r = residu;
}

template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::init( Rep& r) const { return r = zero; }
    
template<class TT> inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::read (Rep& r, const TT a) const { return r = a; }

   // ----- random generators
template<class TT>  template<class RandIter>
inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::nonzerorandom(RandIter& g, Rep& a) const {
    while (iszero( g(a) ) ) {};
    return a;
}

template<class TT>  template<class RandIter>
inline typename OperatorWrapper<TT>::Rep& OperatorWrapper<TT>::random(RandIter& g, Rep& a) const {
    return g(a);
}
