// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id$
// ==========================================================================
// Description:

// ---------
// -- normalized operations
// ---------


// r = a*b
#define __GIVARO_ZPZ64_N_MUL(r,p,a,b) { r = (a*b); r= (r>=p ? r % (uint64)p : r); }
// r *= a
#define __GIVARO_ZPZ64_N_MULIN(r,p,a) { r *= a; r= (r<p ? r : r % (uint64)p);}

// r = a - b
#define __GIVARO_ZPZ64_N_SUB(r,p,a,b) { r = (a-b); r= (r < 0 ? r+p : r);}
// r -= a
#define __GIVARO_ZPZ64_N_SUBIN(r,p,a) { r -= a; r= (r < 0 ? r+p : r);}

// r = a+b
#define __GIVARO_ZPZ64_N_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p);}
// r += a
#define __GIVARO_ZPZ64_N_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p);}

// r <- a*b+c % p
#define __GIVARO_ZPZ64_N_MULADD(r,p,a,b,c) \
{ r = (a*b+c); r= (r<p ? r : r % (uint64)p); }

#define __GIVARO_ZPZ64_N_MULADDIN(r,p,a,b) \
{ r += (a*b); r= (r<p ? r : r % (uint64)p);}

// // a*b-c
// #define __GIVARO_ZPZ64_N_MULSUB(r,p,a,b,c) \
// { r = (a*b+p-c); r= (r<p ? r : r % p); }
// // a*b-c
// #define __GIVARO_ZPZ64_N_MULSUBIN(r,p,a,b) \
// { r -= (a*b)+p; r= (r<p ? r : r % p);}
// a*b-c
#define __GIVARO_ZPZ64_N_MULSUB(r,p,a,b,c) \
{ r = (a*b-c); r= (r<p ? (r>0? r : r% (uint64)p) : r % (uint64)p); }
// a*b-c
#define __GIVARO_ZPZ64_N_MULSUBIN(r,p,a,b) \
{ r -= (a*b); r= (r<p ? (r>0? r : r%(uint64)p) : r % (uint64)p);}


#define __GIVARO_ZPZ64_N_NEG(r,p,a) { r = (a == 0 ? 0 : p-a); }
#define __GIVARO_ZPZ64_N_NEGIN(r,p) { r = (r == 0 ? 0 : p-r); }


inline ZpzDom<Std64>::ZpzDom<Std64>( )
 : _p(0), one(1), zero(0) 
{}

inline ZpzDom<Std64>::ZpzDom<Std64>( Residu_t p )
 : _p(p), one(1), zero(0) 
{}

inline ZpzDom<Std64>::Residu_t ZpzDom<Std64>::residu( ) const
{ return _p; }

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::mul (Rep& r, const Rep a, const Rep b) const
{ 
  register int64 tmp; 
  __GIVARO_ZPZ64_N_MUL(tmp,(int64)_p,(int64)a,(int64)b); 
  return r = (ZpzDom<Std64>::Rep)tmp; 
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::div (Rep& r, const Rep a, const Rep b) const
{ 
  register int64 tmp; 
  register int64 ib; 
  inv(ib, b);
  __GIVARO_ZPZ64_N_MUL(tmp,(int64)_p,(int64)a,(int64)ib); 
  return r = (ZpzDom<Std64>::Rep)tmp; 
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::sub (Rep& r, const Rep a, const Rep b) const
{ 
  register int64 tmp; 
  __GIVARO_ZPZ64_N_SUB(tmp,(int64)_p,(int64)a,(int64)b); 
  return r = (ZpzDom<Std64>::Rep)tmp; 
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::add (Rep& r, const Rep a, const Rep b) const
{ 
  register int64 tmp; 
  __GIVARO_ZPZ64_N_ADD(tmp,(int64)_p,(int64)a,(int64)b); 
  return r = (ZpzDom<Std64>::Rep)tmp; 
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::neg (Rep& r, const Rep a) const
{ 
  register int64 tmp; 
  __GIVARO_ZPZ64_N_NEG(tmp,(int64)_p,(int64)a); 
  return r = (ZpzDom<Std64>::Rep)tmp; 
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::inv (Rep& r, const Rep a) const
{ 
  register int64 d, u, v; 
//   d = ZpzDom<Std64>::gcdext(u, v, a, _p);
  ZpzDom<Std64>::gcdext(d, u, v, a, _p);
//  if ((d != 1) || (d != -1)) GivError::throw_error( GivMathDivZero("Zpz::inv"));
// JGD 04.11.1999
//   if ((d != 1) || (d != -1)) cerr << "GivMathDivZero(Zpz::inv)" << endl;
 // if ((d != 1) && (d != -1)) cerr << "GivMathDivZero(Zpz::inv)" << endl;
  if ((d != 1) && (d != -1)) GivError::throw_error( GivMathDivZero("Zpz::inv"));
// JGD 08.03.2000 : et si u est négatif hein !!!
//   r = (ZpzDom<Std64>::Rep)u; 
//   init(r,u);
  r = (u<0)?(ZpzDom<Std64>::Rep)(u+_p):(ZpzDom<Std64>::Rep)u; 
  if (d == -1) neg(r,r); 
  return r;
}

 // -- inline array operations between ZpzDom<Std64>::Rep
inline void ZpzDom<Std64>::mul (const size_t sz, Array r, constArray a, constArray b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_MUL(tmp, (int64)_p,(int64)a[i], (int64)b[i]); 
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

inline void ZpzDom<Std64>::mul (const size_t sz, Array r, constArray a, Rep b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_MUL(tmp, (int64)_p, (int64)a[i], (int64)b);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

inline void ZpzDom<Std64>::div (const size_t sz, Array r, constArray a, constArray b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    div( r[i], a[i], b[i]); 
  }
}

inline void ZpzDom<Std64>::div (const size_t sz, Array r, constArray a, Rep b) const
{
  ZpzDom<Std64>::Rep ib;
  inv(ib, b);
  mul(sz, r, a, ib);
}

inline void ZpzDom<Std64>::add (const size_t sz, Array r, constArray a, constArray b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_ADD(tmp, (int64)_p, (int64)a[i], (int64)b[i]);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

inline void ZpzDom<Std64>::add (const size_t sz, Array r, constArray a, Rep b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_ADD(tmp, (int64)_p, (int64)a[i], (int64)b);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

inline void ZpzDom<Std64>::sub (const size_t sz, Array r, constArray a, constArray b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_SUB(tmp, (int64)_p, (int64)a[i], (int64)b[i]);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

inline void ZpzDom<Std64>::sub (const size_t sz, Array r, constArray a, Rep b) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_SUB(tmp, (int64)_p, (int64)a[i], (int64)b);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

inline void ZpzDom<Std64>::neg (const size_t sz, Array r, constArray a) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_NEG(tmp, (int64)_p, (int64)a[i]);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}


inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::mulin (Rep& r, const Rep a) const
{ 
  register int64 tmp = (int64)r; 
  __GIVARO_ZPZ64_N_MULIN(tmp,(int64)_p, (int64)a); 
  return r = (ZpzDom<Std64>::Rep)tmp; 
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::divin (Rep& r, const Rep a) const
{ 
  ZpzDom<Std64>::Rep ia;
  inv(ia, a);
  mulin(r, ia);
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::addin (Rep& r, const Rep a) const
{ 
  register int64 tmp = (int64)r; 
  __GIVARO_ZPZ64_N_ADDIN(tmp,(int64)_p, (int64)a); 
  return r = (ZpzDom<Std64>::Rep)tmp; 
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::subin (Rep& r, const Rep a) const
{ 
  register int64 tmp = (int64)r; 
  __GIVARO_ZPZ64_N_SUBIN(tmp,(int64)_p, (int64)a); 
  return r = (ZpzDom<Std64>::Rep)tmp; 
}


inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::negin (Rep& r) const
{ 
  __GIVARO_ZPZ64_N_NEGIN(r,(int64)_p); 
  return r; 
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::invin (Rep& r) const
{ 
  register int64 d, u, v; 
//   d = ZpzDom<Std64>::gcdext(u, v, a, _p);
  ZpzDom<Std64>::gcdext(d, u, v, r, _p);
//  if ((d != 1) || (d != -1)) GivError::throw_error( GivMathDivZero("Zpz::inv"));
// JGD 04.11.1999
//   if ((d != 1) || (d != -1)) cerr << "GivMathDivZero(Zpz::inv)" << endl;
//   if ((d != 1) && (d != -1)) cerr << "GivMathDivZero(Zpz::inv)" << endl;
  if ((d != 1) && (d != -1)) GivError::throw_error( GivMathDivZero("Zpz::invin"));
// JGD 08.03.2000
//   r = (ZpzDom<Std64>::Rep)u; 
  r = (u<0)?(ZpzDom<Std64>::Rep)(u+_p):(ZpzDom<Std64>::Rep)u; 
  if (d == -1) neg(r,r); 
  return r;
}


inline void ZpzDom<Std64>::axpy 
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{ 
  register int64 tmp; 
  __GIVARO_ZPZ64_N_MULADD(tmp, (int64)_p, (int64)a, (int64)b, (int64)c); 
  r = (ZpzDom<Std64>::Rep)tmp; 
}

inline void ZpzDom<Std64>::axpyin 
 (Rep& r, const Rep a, const Rep b) const
{ 
  register int64 tmp = (int64)r; 
  __GIVARO_ZPZ64_N_MULADDIN(tmp, (int64)_p, (int64)a, (int64)b); 
  r = (ZpzDom<Std64>::Rep)tmp; 
}


inline void ZpzDom<Std64>::axpy 
  (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_MULADD(tmp, (int64)_p, (int64)a[i], (int64)x[i], (int64)y[i]);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

inline void ZpzDom<Std64>::axpyin 
  (const size_t sz, Array r, constArray a, constArray x) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp = (int64)r[i];
    __GIVARO_ZPZ64_N_MULADDIN(tmp, (int64)_p, (int64)a[i], (int64)x[i]);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::amxy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
  register int64 tmp;
  __GIVARO_ZPZ64_N_MUL(tmp, (int64)_p, (int64)a, (int64)b);
  __GIVARO_ZPZ64_N_SUB(r, (int64)_p, (int64)c, tmp);
  return r;
}

inline void ZpzDom<Std64>::axmy
 (Rep& r, const Rep a, const Rep b, const Rep c) const
{
  register int64 tmp;
  __GIVARO_ZPZ64_N_MULSUB(tmp, (int64)_p, (int64)a, (int64)b, (int64)c);
  r = (ZpzDom<Std64>::Rep)tmp;
}

// r -= a*b
inline void ZpzDom<Std64>::axmyin 
 (Rep& r, const Rep a, const Rep b) const
{
  register int64 tmp = (int64)r;
  __GIVARO_ZPZ64_N_MULSUBIN(tmp, (int64)_p, (int64)a, (int64)b );
  r = (ZpzDom<Std64>::Rep)tmp;
}


inline void ZpzDom<Std64>::axmy
  (const size_t sz, Array r, constArray a, constArray x, constArray y) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp;
    __GIVARO_ZPZ64_N_MULSUB(tmp, (int64)_p, (int64)a[i], (int64)x[i], (int64)y[i]);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}

// r -= a*b
inline void ZpzDom<Std64>::axmyin
  (const size_t sz, Array r, constArray a, constArray x) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    register int64 tmp = (int64)r[i];
    __GIVARO_ZPZ64_N_MULSUBIN(tmp, (int64)_p, (int64)a[i], (int64)x[i]);
    r[i] = (ZpzDom<Std64>::Rep)tmp;
  }
}


// ---------
// -- misc operations
// ---------

inline void ZpzDom<Std64>::assign 
  ( const size_t sz, Array r, constArray a ) const
{
  for (register size_t i=sz-1; i!=0; --i) {
    if (a[i] <ZpzDom<Std64>::zero) {
       r[i] = a[i] + _p;
       if (r[i] <ZpzDom<Std64>::zero) r[i] = r[i] % (uint64)_p;
    }
    else if (a[i] >_p) {
       r[i] = a[i] - _p;
       if (r[i] >=_p) r[i] = r[i] % (uint64)_p;
    }
    else r[i] = a[i];
  }
}

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign ( Rep& r, const long a ) const
{  return init(r, a);
}

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const unsigned long a ) const
{ return r = (Rep)( a >= (unsigned long)_p ? a % (unsigned long)_p : a); 
}


inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const long a ) const
{
  int sign; unsigned long ua;
  if (a <0) { sign =-1; ua = -a;}
  else { ua = a; sign =1; }
  r = (ua >=_p) ? ua % (uint64)_p : ua;
  if (r && (sign ==-1)) r = _p - r;
  return r;
}

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const unsigned long long a ) const
{ return r = (Rep)( a >= (unsigned long)_p ? a % (unsigned long)_p : a); 
}


inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::init ( Rep& r, const long long a ) const
{
  int sign; unsigned long ua;
  if (a <0) { sign =-1; ua = -a;}
  else { ua = a; sign =1; }
  r = (ua >=_p) ? ua % (uint64)_p : ua;
  if (r && (sign ==-1)) r = _p - r;
  return r;
}

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign ( Rep& r, const int a ) const
{ return ZpzDom<Std64>::assign( r, (long)a); }

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign ( Rep& r, const unsigned long a ) const
{ return init(r,a); }

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign 
  ( Rep& r, const unsigned int a ) const
{ return init(r, (unsigned long)a); }

inline  ZpzDom<Std64>::Rep&  ZpzDom<Std64>::assign 
  ( Rep& r, const Rep a ) const
{ return assign(r, (long)a); }

// inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::GIVRANDOM(Rep& a) const {
// 	assign(a, lrand48());
// }

// inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::NONZEROGIVRANDOM(Rep& a) const {
// 	while (iszero(assign(a, lrand48()))) {};
// }

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::random(RandIter& g, Rep& a) const {
	init(a, g());
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::random(RandIter& g, Rep& a, const Rep& b) const {
	init(a, g());
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::random(RandIter& g, Rep& a, long b) const {
	init(a, g() %(uint64) b);
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::nonzerorandom(RandIter& g, Rep& a) const {
	while (iszero(init(a, g()))) {};
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::nonzerorandom(RandIter& g, Rep& a, const Rep& b) const {
	while (iszero(init(a, g()))) {};
}

template< class RandIter >
inline  ZpzDom<Std64>::Rep& ZpzDom<Std64>::nonzerorandom(RandIter& g, Rep& a, long b) const {
	while (iszero(init(a, g() %(uint64) b))) {};
}

inline void ZpzDom<Std64>::init 
  ( const size_t sz, Array r, constArray a ) const
{
  for (register size_t i=sz-1; i!=0; --i) 
       r[i] = a[i];
}

inline ZpzDom<Std64>::Rep& ZpzDom<Std64>::init ( Rep& r ) const
{ return r = zero; }

inline void ZpzDom<Std64>::dotprod 
  ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const
{
  int stride = 1;
  if ((unsigned long)bound < GIVARO_MAXUINT64)
   stride = GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound);
  unsigned long dot = zero;
  if ((sz <10) && (sz <stride)) {
    for( register int i= sz-1; i>=0; --i) 
      dot += a[i] * b[i]; 
    if (dot > _p) r = (Rep)(dot % (uint64)_p);
    else r = (Rep)dot;
    return;
  }
  int i_begin=0;
  stride &= ~0x1;
  if (stride ==0) {
    for( register int i= sz-1; i>0; --i) {
      dot += a[i] * b[i];
      if (dot>_p) dot %= _p;
    }
    r = (Rep)dot;
    return;
  }
  do {
    size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
    if (min_sz & 0x1 !=0) 
      { min_sz--; i_begin++; dot += a++[min_sz] * b++[min_sz]; }
    if (min_sz > 1) 
      for( register size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b ) 
      {
        dot += a[0] * b[0]; 
        dot += a[1] * b[1];
      }
    if (dot>(int64)_p) dot %= (uint64)_p;
    i_begin += min_sz;
  } while (i_begin <sz);
  r = (Rep)dot;
}

inline void ZpzDom<Std64>::dotprod
  ( Rep& r, const size_t sz, constArray a, constArray b ) const
{
  return ZpzDom<Std64>::dotprod(r, _p, sz, a, b);
}


  //  a -> r: int64 to double
inline void 
  ZpzDom<Std64>::i2d ( const size_t sz, double* r, constArray a ) const
{
  for (size_t i=0; i<sz; ++i) r[i] = a[i];
}

  //  a -> r: double to int64 
inline void 
  ZpzDom<Std64>::d2i ( const size_t sz, Array r, const double* a ) const
{
  union d_2_l {
    double d;
    int64 r[2];
  };
//  static const double offset = 4503599627370496.0; // 2^52
  double offset = 4503599627370496.0; // 2^52
  for (size_t i=0; i<sz; ++i)
  {
      register d_2_l tmp;
      // - normalization: put fractional part at the end of the representation
      tmp.d = a[i] + offset; 
      r[i] = tmp.r[1];
      if (r[i] <_p) r[i] %= _p;
  }
  //    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]-_p);
  //    r[i] = (r[i] <_p ? r[i] : r[i]%_p);
  //    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]%_p);
}


 // ------------------------- Miscellaneous functions

inline int ZpzDom<Std64>::areEqual(const Rep a, const Rep b) const
{ return a == b; }

inline int ZpzDom<Std64>::areNEqual(const Rep a, const Rep b) const
{ return a != b; }

inline int ZpzDom<Std64>::iszero(const Rep a) const
{ return a == ZpzDom<Std64>::zero; }

inline int ZpzDom<Std64>::isnzero(const Rep a) const
{ return a != ZpzDom<Std64>::zero; }

inline int ZpzDom<Std64>::isone(const Rep a) const
{ return a == ZpzDom<Std64>::one; }

inline size_t ZpzDom<Std64>::length(const Rep a) const
{ return ZpzDom<Std64>::size_rep;}

 // -- Input: (z, <_p>)
inline istream& ZpzDom<Std64>::read (istream& s) 
{
  char ch; 
  s >> ws >> ch;
  if (ch != '(')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std64>::read: syntax error: no '('"));
    cerr << "GivBadFormat(ZpzDom<Std64>::read: syntax error: no '('))" << endl;

  s >> ws >> ch;
  if (ch != 'z')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std64>::read: bad domain object"));
    cerr << "GivBadFormat(ZpzDom<Std64>::read: bad domain object))" << endl;

  s >> ws >> ch;
  if (ch != ',')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std64>::read: syntax error: no ','"));
    cerr << "GivBadFormat(ZpzDom<Std64>::read: syntax error: no ',')) " << endl;

  s >> ws >> _p;

  s >> ws >> ch;
  if (ch != ')')
//    GivError::throw_error( GivBadFormat("ZpzDom<Std64>::read: syntax error: no ')'"));
    cerr << "GivBadFormat(ZpzDom<Std64>::read: syntax error: no ')')) " << endl;

  return s;
}

inline ostream& ZpzDom<Std64>::write (ostream& s ) const
{
  return s << "(z," << residu() << ')';
}

inline istream& ZpzDom<Std64>::read (istream& s, Rep& a) const
{
  s >> a;
  init(a, a);
  return s;
}

inline ostream& ZpzDom<Std64>::write (ostream& s, const Rep a) const
{
  return s << a;
}


inline Integer& ZpzDom<Std64>::write (Integer& r, const Rep a) const
{
  return r = Integer(a);
}
