// ==========================================================================
// file: lin_dom_tt.inl
// Description:
// Bugs:
// Time-stamp: <05 Apr 00 13:13:44 Jean-Guillaume.Dumas@imag.fr> 
// ==========================================================================

 // ------------------------- Miscellaneous functions

template<class TT> inline bool TTDom<TT>::iszero(const Rep& a) const
  { return a == zero ; }

template<class TT> inline bool TTDom<TT>::isone(const Rep& a) const
  { return a == one ; }

template<class TT> inline bool TTDom<TT>::isnzero(const Rep& a) const
  { return a != zero ; }

template<class TT> inline bool TTDom<TT>::isnone(const Rep& a) const
  { return a != one ; }

template<class TT> inline size_t TTDom<TT>::length(const Rep& ) const
  { return sizeof(TT) ;}

template<class TT> inline     bool TTDom<TT>::isequal( const Rep& a, const Rep& b) const { return a == b ; }
    
template<class TT> inline     bool TTDom<TT>::isnequal( const Rep& a, const Rep& b) const { return a != b ; }
    
template<class TT> inline     bool TTDom<TT>::islt( const Rep& a, const Rep& b) const { return a < b ; }
    
template<class TT> inline     bool TTDom<TT>::isgt( const Rep& a, const Rep& b) const { return a > b ; }
    
  // ----------- Usefull method :
template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::mul (Rep& r, const Rep&a, const Rep& b) const 
  { return r = a * b; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::mulin (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a) const 
  { return r *= a; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::div (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a, const TTDom<TT>::Rep&b) const 
  { return r = a / b; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::divin (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a) const 
  { return r /= a; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::add (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a, const TTDom<TT>::Rep&b) const 
  { return r = a + b; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::addin (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a) const 
  { return r += a; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::sub (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a, const TTDom<TT>::Rep&b) const 
  { return r = a - b; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::subin (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a) const 
  { return r -= a; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::neg (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a) const 
  { return r = -a; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::negin (TTDom<TT>::Rep& r) const 
  { return r = -r; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::inv (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a) const 
  { return r = one/a; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::invin (TTDom<TT>::Rep& r) const 
  { return r = one/r; }

/*
#include "lin_sqrt.h"

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::sqrt (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a) const 
  { return r = ::sqrt(a); }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::sqrtin (TTDom<TT>::Rep& r) const 
  { return r = ::sqrt(r); }
*/

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::axpy (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a, const TTDom<TT>::Rep&b, const TTDom<TT>::Rep&c) const 
  { return r = a * b + c; }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::axpyin (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a, const TTDom<TT>::Rep&b) const { return r = r + a*b ; }

  // -- axmyin: r <- r - a * b mod p
template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::axmyin (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a, const TTDom<TT>::Rep&b) const { return r -= a*b;  }

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::axmy (TTDom<TT>::Rep& r, const TTDom<TT>::Rep&a, const TTDom<TT>::Rep&b, const TTDom<TT>::Rep&c)  const { return r = a * b - c; }


template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::assign( Rep& r, const Rep&residu ) const {
    return r = residu;
}

template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::init( Rep& r) const { return r = zero; }
    
template<class TT> inline TTDom<TT>::Rep& TTDom<TT>::read (Rep& r, const TT a) const { return r = a; }

   // ----- random generators
template<class TT>  template<class RandIter>
inline TTDom<TT>::Rep& TTDom<TT>::nonzerorandom(RandIter& g, Rep& a) const {
    while (iszero( g(a) ) ) {};
    return a;
}

template<class TT>  template<class RandIter>
inline TTDom<TT>::Rep& TTDom<TT>::random(RandIter& g, Rep& a) const {
    return g(a);
}
