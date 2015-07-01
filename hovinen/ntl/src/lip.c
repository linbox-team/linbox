
#include <NTL/lip.h>
#include <NTL/IsFinite.h>


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if (defined(NTL_CPLUSPLUS_ONLY) && !defined(__cplusplus))
#error "CPLUSPLUS_ONLY flag set...must use C++ compiler"
#endif


#define DENORMALIZE  NTL_FDOUBLE_PRECISION


#define MIN_SETL	(4)
   /* zsetlength allocates a multiple of MIN_SETL digits */

#define MustAlloc(c, len)  (!(c) || ((c)[-1] >> 1) < (len))
   /* Fast test to determine if allocation is necessary */


#if (defined(NTL_SINGLE_MUL) && !defined(NTL_FAST_INT_MUL))
#define MulLo(rres,a,b) \
{ \
   double _x = ((double) a) * ((double) b) + (DENORMALIZE); \
   unsigned long _x_lo; \
   NTL_FetchLo(_x_lo, _x); \
   rres =  _x_lo; \
}
#else
#define MulLo(rres,a,b) rres = (a)*(b)
#endif

#if (defined(NTL_SINGLE_MUL))

#define zaddmulp(a,b,d,t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   unsigned long  __lhi, __llo;\
   double __lx = ((double) ((_a) + (_t))) + ((double) _b)*((double) _d); \
  __lx += (DENORMALIZE); \
   NTL_FetchHiLo(__lhi,__llo,__lx);\
   __lhi = ((__lhi<<6)|(__llo>>26)) & 0x3FFFFFF; \
   __llo &= 0x3FFFFFF; \
   (a) = __llo;\
   (t) = __lhi;\
}

#define zxmulp(a,b,d,t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   unsigned long  __lhi, __llo;\
   double __lx = ((double) ((_t))) + ((double) _b)*((double) _d); \
  __lx += (DENORMALIZE); \
   NTL_FetchHiLo(__lhi,__llo,__lx);\
   __lhi = ((__lhi<<6)|(__llo>>26)) & 0x3FFFFFF; \
   __llo &= 0x3FFFFFF; \
   (a) = __llo;\
   (t) = __lhi;\
}

#define zaddmulpsq(a,b,t) \
{ \
   long _a = (a), _b = (b); \
   unsigned long  __lhi, __llo; \
   double __lx = ((double) (_a)) + ((double) _b)*((double) _b); \
   __lx += (DENORMALIZE); \
   NTL_FetchHiLo(__lhi,__llo,__lx);\
   __lhi = ((__lhi<<6)|(__llo>>26)) & 0x3FFFFFF; \
   __llo &= 0x3FFFFFF; \
   (a) =  __llo;\
   (t) =  __lhi;\
}

#else

#if (defined(NTL_LONG_LONG))

#ifndef NTL_LONG_LONG_TYPE
#define NTL_LONG_LONG_TYPE long long
#endif

#define zaddmulp(a, b, d, t) { \
   NTL_LONG_LONG_TYPE pp = ((NTL_LONG_LONG_TYPE) (b)) * ((NTL_LONG_LONG_TYPE) (d)) + ((t)+(a)); \
   (a) = ((long)(pp)) & NTL_RADIXM; \
   (t) = (long) (pp >> NTL_NBITS); \
} 


#define zxmulp(a, b, d, t) { \
   NTL_LONG_LONG_TYPE pp = ((NTL_LONG_LONG_TYPE) (b)) * ((NTL_LONG_LONG_TYPE) (d)) + (t); \
   (a) = ((long)(pp)) & NTL_RADIXM; \
   (t) = (long) (pp >> NTL_NBITS); \
} 

#define zaddmulpsq(a,b,t) { \
   NTL_LONG_LONG_TYPE pp = ((NTL_LONG_LONG_TYPE) (b)) * ((NTL_LONG_LONG_TYPE) (b)) + (a); \
   (a) = ((long)(pp)) & NTL_RADIXM; \
   (t) = (long) (pp >> NTL_NBITS); \
}


#elif (defined(NTL_AVOID_FLOAT))
   

#define zaddmulp(  a,  b,  d,  t) { \
        unsigned long b1 = b & NTL_RADIXROOTM; \
        unsigned long d1 = d & NTL_RADIXROOTM; \
        unsigned long bd,b1d1,m,aa= (a) + (t); \
	unsigned long ld = (d>>NTL_NBITSH); \
	unsigned long lb = (b>>NTL_NBITSH); \
 \
        bd=lb*ld; \
        b1d1=b1*d1; \
        m=(lb+b1)*(ld+d1) - bd - b1d1; \
        aa += ( b1d1+ ((m&NTL_RADIXROOTM)<<NTL_NBITSH)); \
        (t) = (aa >> NTL_NBITS) + bd + (m>>NTL_NBITSH); \
        (a) = aa & NTL_RADIXM; \
}





#define zxmulp(  a,  b,  d,  t) { \
        unsigned long b1 = b & NTL_RADIXROOTM; \
        unsigned long d1 = d & NTL_RADIXROOTM; \
        unsigned long bd,b1d1,m,aa= (t); \
	unsigned long ld = (d>>NTL_NBITSH); \
	unsigned long lb = (b>>NTL_NBITSH); \
 \
        bd=lb*ld; \
        b1d1=b1*d1; \
        m=(lb+b1)*(ld+d1) - bd - b1d1; \
        aa += ( b1d1+ ((m&NTL_RADIXROOTM)<<NTL_NBITSH)); \
        (t) = (aa >> NTL_NBITS) + bd + (m>>NTL_NBITSH); \
        (a) = aa & NTL_RADIXM; \
}




#define zaddmulpsq(_a, _b, _t) \
{ \
	long lb = (_b); \
	long b1 = (_b) & NTL_RADIXROOTM; \
	long aa = (_a) + b1 * b1; \
 \
	b1 = (b1 * (lb >>= NTL_NBITSH) << 1) + (aa >> NTL_NBITSH); \
	aa = (aa & NTL_RADIXROOTM) + ((b1 & NTL_RADIXROOTM) << NTL_NBITSH); \
	(_t) = lb * lb + (b1 >> NTL_NBITSH) + (aa >> NTL_NBITS); \
	(_a) = (aa & NTL_RADIXM); \
}



#else

#if (NTL_BITS_PER_LONG <= NTL_NBITS + 2)

#if (NTL_ARITH_RIGHT_SHIFT)
/* value right-shifted is -1..1 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ); \
   _t2 = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   _t1 = (_t1 & NTL_RADIXM) + _a +_t; \
   (t) = _t2 + (((unsigned long)_t1) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}


#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/* value shifted is -1..1 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b); \
   long _t1 = _b*_b; \
   long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ); \
   _t2 = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   _t1 = (_t1 & NTL_RADIXM) + _a; \
   (t) = _t2 + (((unsigned long)_t1) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zam_decl double _ds; long _hi, _lo, _s;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = ((_s << 1)+1)*(NTL_FRADIX_INV/2.0); \
   _lo = _b*_s; \
   _hi = (long) (((double) _b)*_ds); \
}

/* value shifted is 0..3 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is -1..+1 */
#define zsx_loop(a,t,nb) \
{ \
   long  _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}


/* value shifted is -2..+1 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _a + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..3 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* value shifted is -1..+1 */
#define zsx_finish(a,t) \
{ \
   long _t = (t); \
   _lo = _lo + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* value shifted is -2..+1 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _a + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
}

#else

/* right shift is not arithmetic */

/* value right-shifted is 0..2 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   _t2 = _t2 + ( ((unsigned long) (_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS ); \
   _t1 = (_t1 & NTL_RADIXM) + _a +_t; \
   (t) = _t2 + (((unsigned long)_t1) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}


#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/* value shifted is 0..2 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b); \
   long _t1 = _b*_b; \
   long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ) - 1; \
   _t2 = _t2 + ( ((unsigned long) (_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS ); \
   _t1 = (_t1 & NTL_RADIXM) + _a; \
   (t) = _t2 + (((unsigned long)_t1) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zam_decl double _ds; long _hi, _lo, _s;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = ((_s << 1)+1)*(NTL_FRADIX_INV/2.0); \
   _lo = _b*_s; \
   _hi = (long) (((double) _b)*_ds); \
}

/* value shifted is 0..3 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..2 */
#define zsx_loop(a,t,nb) \
{ \
   long _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..3 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _hi += 2; \
   _lo = _a + _t - _lo; \
   (t) = (((unsigned long)(_lo + (_hi<<NTL_NBITS))) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..3 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* value shifted is 0..2 */
#define zsx_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* value shifted is 0..3 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _hi += 2; \
   _lo = _a + _t - _lo; \
   (t) = (((unsigned long)(_lo + (_hi<<NTL_NBITS))) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
}

#endif
   
#else

/*  NTL_BITS_PER_LONG > NTL_NBITS + 2, and certain optimizations can be
    made.  Useful on 64-bit machines.  */

#if (NTL_ARITH_RIGHT_SHIFT)
/* shift is -1..+3 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _a + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ); \
   (t) = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ); \
   (t) = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/* shift is -1..+2 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b), _t = (t); \
   long _t1 = _b*_b + _a; \
   long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ); \
   (t) = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zam_decl double _ds; long _hi, _lo, _s;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = _s*NTL_FRADIX_INV; \
   _lo = _b*_s; \
   _hi = (long) (((double) _b)*_ds); \
}

/* shift is -1..+3 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _a + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is -1..+2 */
#define zsx_loop(a,t,nb) \
{ \
   long _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is -3..+1 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _a + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is -1..+3 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _a + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is -1..+2 */
#define zsx_finish(a,t) \
{ \
   long _t = (t); \
   _lo = _lo + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is -3..+1 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _a + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
}

#else

/* right shift is not arithmetic */

/* shift is 0..4 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _a + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/* shift is 0..3 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b), _t = (t); \
   long _t1 = _b*_b + _a; \
   long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zam_decl double _ds; long _hi, _lo, _s;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = _s*NTL_FRADIX_INV; \
   _lo = _b*_s; \
   _hi = (long) (((double) _b)*_ds); \
}

/* shift is 0..4 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _hi--; \
   _lo = _lo + _a + _t; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..3 */
#define zsx_loop(a,t,nb) \
{ \
   long _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _hi--; \
   _lo = _lo + _t; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..4 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _hi += 3; \
   _lo = _a + _t - _lo; \
   (t) = (((unsigned long)(_lo + (_hi<<NTL_NBITS))) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..4 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is 0..3 */
#define zsx_finish(a,t) \
{ \
   long _t = (t); \
   _lo = _lo + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is 0..4 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _hi += 3; \
   _lo = _a + _t - _lo; \
   (t) = (((unsigned long)(_lo + (_hi<<NTL_NBITS))) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
}
#endif

#endif

#endif

#endif

static
void zaddmulone(long *lama, long *lamb)
{ 
   long lami; 
   long lams = 0; 
 
   lams = 0; 
   for (lami = (*lamb++); lami > 0; lami--) { 
      lams += (*lama + *lamb++); 
      *lama++ = lams & NTL_RADIXM; 
      lams >>= NTL_NBITS; 
   } 
   *lama += lams; 
}

#if (NTL_ARITH_RIGHT_SHIFT)
static
void zsubmulone(long *lama, long *lamb)
{ 
   long lami; 
   long lams = 0; 
 
   lams = 0; 
   for (lami = (*lamb++); lami > 0; lami--) { 
      lams += (*lama - *lamb++); 
      *lama++ = lams & NTL_RADIXM; 
      lams >>= NTL_NBITS; 
   } 
   *lama += lams; 
}
#else
static
void zsubmulone(long *lama, long *lamb)
{ 
   long lami; 
   long lams = 0; 
 
   lams = 0; 
   for (lami = (*lamb++); lami > 0; lami--) { 
      lams = *lama - *lamb++ - lams; 
      *lama++ = lams & NTL_RADIXM; 
      lams = (lams < 0);
   } 
   *lama -= lams; 
}

#endif

#if (defined(NTL_SINGLE_MUL))

void zaddmul(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = (*amb++);
   double dams = (double) ams;
   double xx;
   double yy;
   unsigned long lo_wd, lo;
   unsigned long hi_wd, hi;

   xx  =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams +DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = lo + (*ama) + carry;
      *ama = lo & 0x3FFFFFF;
      carry = hi + (lo >> 26);
      ama++; 
      xx = yy;
   } 

   NTL_FetchHiLo(hi_wd, lo_wd, xx);
   lo = lo_wd & 0x3FFFFFF;
   hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
   lo = lo + (*ama) + carry;
   *ama = lo & 0x3FFFFFF;
   carry = hi + (lo >> 26);
   ama++; 
   *ama += carry; 
}

static
void zsxmul(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = (*amb++);
   double dams = (double) ams;
   double xx;
   double yy;
   unsigned long lo_wd, lo;
   unsigned long hi_wd, hi;

   xx  =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams +DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = lo + carry;
      *ama = lo & 0x3FFFFFF;
      carry = hi + (lo >> 26);
      ama++; 
      xx = yy;
   } 

   NTL_FetchHiLo(hi_wd, lo_wd, xx);
   lo = lo_wd & 0x3FFFFFF;
   hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
   lo = lo + carry;
   *ama = lo & 0x3FFFFFF;
   carry = hi + (lo >> 26);
   ama++; 
   *ama = carry; 
}

static
void zaddmulsq(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = ams;
   double dams = (double) (*amb++);
   double xx;
   double yy;
   unsigned long lo, lo_wd;
   unsigned long hi, hi_wd;

   xx =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams + DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = lo + (*ama) + carry;
      *ama = lo & 0x3FFFFFF;
      carry = hi + (lo >> 26);
      ama++; 
      xx = yy;
   } 
   if (i==1) {
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = lo + (*ama) + carry;
      *ama = lo & 0x3FFFFFF;
      carry = hi + (lo >> 26);
      ama++; 
   }
   *ama += carry; 
}

#else

#if (defined(NTL_AVOID_FLOAT) || defined(NTL_LONG_LONG))
void zaddmul(long lams, long *lama, long *lamb)
{
        long lami;
        long lamcarry = 0;

        for (lami = (*lamb++); lami > 0; lami--)
        {
                zaddmulp(*lama, *lamb, lams, lamcarry);
                lama++;
                lamb++;
        }
        *lama += lamcarry;
}

static
void zsxmul(long lams, long *lama, long *lamb)
{
        long lami;
        long lamcarry = 0;

        for (lami = (*lamb++); lami > 0; lami--)
        {
                zxmulp(*lama, *lamb, lams, lamcarry);
                lama++;
                lamb++;
        }
        *lama = lamcarry;
}

static
void zaddmulsq(long lsqi, long *lsqa, long *lsqb)
{
        long lsqs = *(lsqb);
        long lsqcarry = 0;

        lsqb++;
        for (; lsqi > 0; lsqi--)
        {
                zaddmulp(*lsqa, *lsqb, lsqs, lsqcarry);
                lsqa++;
                lsqb++;
        }
        *lsqa += lsqcarry;
}



#else

void zaddmul(long lams, long *lama, long *lamb)
{ 
   long lami = (*lamb++)-1; 
   long lamcarry = 0; 
   zam_decl;

   zam_init(*lamb, lams);
   lamb++;

 
   for (; lami > 0; lami--) { 
      zam_loop(*lama, lamcarry, *lamb);
      lama++; 
      lamb++; 
   } 
   zam_finish(*lama, lamcarry);
   lama++;
        *lama += lamcarry; 
}



static
void zsxmul(long lams, long *lama, long *lamb)
{ 
   long lami = (*lamb++)-1; 
   long lamcarry = 0; 
   zam_decl;

   zam_init(*lamb, lams);
   lamb++;

 
   for (; lami > 0; lami--) { 
      zsx_loop(*lama, lamcarry, *lamb);
      lama++; 
      lamb++; 
   } 
   zsx_finish(*lama, lamcarry);
   lama++;
        *lama = lamcarry; 
}



static
void zaddmulsq(long lsqi, long *lsqa, long *lsqb)
{ 
   long lsqs; 
   long lsqcarry; 
   zam_decl

   if (lsqi <= 0) return;

   lsqs = *lsqb;
   lsqcarry = 0;

   lsqb++; 
   zam_init(*lsqb, lsqs);
   lsqb++;
   lsqi--;
   for (; lsqi > 0; lsqi--) { 
      zam_loop(*lsqa, lsqcarry, *lsqb);
      lsqa++; 
      lsqb++; 
   } 
   zam_finish(*lsqa, lsqcarry);
   lsqa++;
   *lsqa += lsqcarry; 
}


#endif
#endif


#if (defined(NTL_SINGLE_MUL))

#if (NTL_ARITH_RIGHT_SHIFT)

static
void zsubmul(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = (*amb++);
   double dams = (double) ams;
   double xx;
   double yy;
   unsigned long lo_wd, lo;
   unsigned long hi_wd, hi;

   xx  =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams +DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = (*ama) + carry - lo;
      *ama = lo & 0x3FFFFFF;
      carry = (((long)lo) >> 26) - hi;
      ama++; 
      xx = yy;
   } 

   NTL_FetchHiLo(hi_wd, lo_wd, xx);
   lo = lo_wd & 0x3FFFFFF;
   hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
   lo = (*ama) + carry - lo;
   *ama = lo & 0x3FFFFFF;
   carry = (((long)lo) >> 26) - hi;
   ama++; 
   *ama += carry; 
}

#else

static
void zsubmul(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = (*amb++);
   double dams = (double) ams;
   double xx;
   double yy;
   unsigned long lo_wd, lo;
   unsigned long hi_wd, hi;

   xx  =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams +DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = (*ama) + carry - lo;
      *ama = lo & 0x3FFFFFF;
      carry = ((lo + (1L << 27)) >> 26) - hi - 2;
      ama++; 
      xx = yy;
   } 

   NTL_FetchHiLo(hi_wd, lo_wd, xx);
   lo = lo_wd & 0x3FFFFFF;
   hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
   lo = (*ama) + carry - lo;
   *ama = lo & 0x3FFFFFF;
   carry = ((lo + (1L << 27)) >> 26) - hi - 2;
   ama++; 
   *ama += carry; 
}

#endif

#else

#if (defined(NTL_LONG_LONG))

static
void zsubmul(long lams, long *lama, long *lamb)
{
        long lami;
        long lamcarry = 0;

        lams = -lams;

        for (lami = (*lamb++); lami > 0; lami--)
        {
                zaddmulp(*lama, *lamb, lams, lamcarry);
                lama++;
                lamb++;
        }
        *lama += lamcarry;
}

#elif (defined(NTL_AVOID_FLOAT))

static void
zsubmul(
        long r,
        verylong a,
        verylong b
        )
{
        long rd = NTL_RADIX - r;
        long i;
        long carry = NTL_RADIX;

        for (i = (*b++); i > 0; i--)
        {
                zaddmulp(*a, *b, rd, carry);
                a++;
                carry += NTL_RADIXM - (*b++);
        }
        *a += carry - NTL_RADIX; /* unnormalized */
}



#else

static
void zsubmul(long lams, long *lama, long *lamb)
{ 
   long lami = (*lamb++)-1; 
   long lamcarry = 0; 
   zam_decl;

   zam_init(*lamb, lams);
   lamb++;

 
   for (; lami > 0; lami--) { 
      zam_subloop(*lama, lamcarry, *lamb);
      lama++; 
      lamb++; 
   } 
   zam_subfinish(*lama, lamcarry);
   lama++;
   *lama += lamcarry; 
}

#endif

#endif




/*
	zdiv21 returns quot, numhigh so

	quot = (numhigh*NTL_RADIX + numlow)/denom;
	numhigh  = (numhigh*NTL_RADIX + numlow)%denom;

Assumes 0 <= numhigh < denom < NTL_RADIX and 0 <= numlow < NTL_RADIX.
*/



#define zdiv21(numhigh, numlow, denom, deninv, quot) \
{ \
   long lr21; \
   long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) \
          + (double) (numlow)) * (deninv)); \
   long lp21; \
   MulLo(lp21, lq21, denom); \
   lr21 = (numhigh << NTL_NBITS) + numlow - lp21; \
   if (lr21 < 0) { \
      lq21--; \
      lr21 += denom; \
   } \
   else if (lr21 >= denom) { \
      lr21 -= denom; \
      lq21++; \
   } \
   quot = lq21; \
   numhigh = lr21; \
}


#if 0

#define zdiv21(numhigh, numlow, denom, deninv, quot) \
{ \
   long lr21; \
   long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) \
          + (double) (numlow)) * (deninv)); \
   long lp21; \
   MulLo(lp21, lq21, denom); \
   lr21 = (numhigh << NTL_NBITS) + numlow - lp21; \
   if (lr21 < 0) \
   { \
      do \
      { \
         lq21--; \
      } while ((lr21 += denom) < 0); \
   } \
   else \
   { \
      while (lr21 >= denom) \
      { \
         lr21 -= denom; \
         lq21++; \
      }; \
   } \
   quot = lq21; \
   numhigh = lr21; \
}

#endif

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
#define zrem21(numhigh, numlow, denom, deninv) \
{ \
   long lr21; \
   long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) \
          + (double) (numlow)) * (deninv)); \
   long lp21; \
   MulLo(lp21, lq21, denom); \
   lr21 = (numhigh << NTL_NBITS) + numlow - lp21; \
   lr21 += (lr21 >> (NTL_BITS_PER_LONG-1)) & denom; \
   lr21 -= denom; \
   lr21 += (lr21 >> (NTL_BITS_PER_LONG-1)) & denom; \
   numhigh = lr21; \
}
#else
#define zrem21(numhigh, numlow, denom, deninv) \
{ \
   long lr21; \
   long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) \
      + (double) (numlow)) * (deninv)); \
   long lp21; \
   MulLo(lp21, lq21, denom); \
   lr21 = (numhigh << NTL_NBITS) + numlow - lp21; \
   if (lr21 < 0) lr21 += denom; \
   else if (lr21 >= denom) lr21 -= denom; \
   numhigh = lr21; \
}
#endif



static 
void zhalt(char *c)
{
   fprintf(stderr,"fatal error:\n   %s\nexit...\n",c);
   fflush(stderr);
   abort();
}


void zsetlength(verylong *v, long len)
{
   verylong x = *v;

   if (len < 0)
      zhalt("negative size allocation in zsetlength");

   if (len >= (1L << (NTL_BITS_PER_LONG-4))/NTL_NBITS)
      zhalt("size too big in zsetlength");

   if (x) {
      long oldlen = x[-1];
      long fixed = oldlen & 1;

      oldlen = oldlen >> 1;

      if (fixed) {
         if (len > oldlen) 
            zhalt("internal error: can't grow this verylong");
         else
            return;
      }

      if (len <= oldlen) return;

      len++;  /* always allocate at least one more than requested */

      oldlen = (long) (oldlen * 1.2); /* always increase by at least 20% */
      if (len < oldlen)
         len = oldlen;

      /* round up to multiple of MIN_SETL */
      len = ((len+(MIN_SETL-1))/MIN_SETL)*MIN_SETL; 

      x[-1] = len << 1;
      if (!(x = (verylong)realloc((void*)(&(x[-1])), (len + 2) * (sizeof (long))))) {
         zhalt("reallocation failed in zsetlength");
      }
   }
   else {
      len++;
      len = ((len+(MIN_SETL-1))/MIN_SETL)*MIN_SETL;

      if (!(x = (verylong)malloc((len + 2)*(sizeof (long))))) {
         zhalt("allocation failed in zsetlength");
      }
      x[0] = len << 1;
      x[1] = 1;
      x[2] = 0;
   }

   *v = x+1;
}

void zfree(verylong *x)
{
   if (!(*x))
      return;

   if ((*x)[-1] & 1)
      zhalt("Internal error: can't free this verylong");

   {
      verylong y = (*x - 1);
      free((void*)y);
      *x = 0;
      return;
   }
}

double zdoub(verylong n)
{
   double res;
   long i;

   if (!n)
      return ((double) 0);
   if ((i = n[0]) < 0)
      i = -i;
   res = (double) (n[i--]);
   for (; i; i--)
      res = res * NTL_FRADIX + (double) (n[i]);
   if (n[0] > 0)
      return (res);
   return (-res);
}



void zdoubtoz(double a, verylong *xx)
{
   verylong x;
   long neg, i, t, sz;


   a = floor(a);

   if (!IsFinite(&a))
      zhalt("zdoubtoz: attempt to convert non-finite value");


   if (a < 0) {
      a = -a;
      neg = 1;
   }
   else
      neg = 0;

   if (a == 0) {
      zzero(xx);
      return;
   }

   sz = 1;
   a = a*NTL_FRADIX_INV;

   while (a >= 1) {
      a = a*NTL_FRADIX_INV;
      sz++;
   }
         
   x = *xx;
   if (MustAlloc(x, sz)) {
      zsetlength(&x, sz);
      *xx = x;
   }

   for (i = sz; i > 0; i--) {
      a = a*NTL_FRADIX;
      t = (long) a;
      x[i] = t;
      a = a - t;
   }

   x[0] = (neg ? -sz : sz);
}

void zzero(verylong *aa)
{
   if (!(*aa)) zsetlength(aa, 1);
   (*aa)[0] =  1;
   (*aa)[1] =  0;
}

void zone(verylong *aa)
{
   if (!(*aa)) zsetlength(aa, 1);
   (*aa)[0] =  1;
   (*aa)[1] =  1;
}

void zcopy(verylong a, verylong *bb)
{
   long i;
   verylong b = *bb;

   if (!a) {
      zzero(bb);
      return;
   }
   if (a != b) {
      if ((i = *a) < 0)
         i = (-i);
      if (MustAlloc(b, i)) {
         zsetlength(&b, i);
         *bb = b;
      }
      for (; i >= 0; i--)
         *b++ = *a++;
   }
}

void zintoz(long d, verylong *aa)
{
   long i;
   long anegative;
   unsigned long d1, d2;
   verylong a = *aa;

   anegative = 0;
   if (d < 0) {
      anegative = 1;
      d1 = -d;
   }
   else
      d1 = d;

   i = 0;
   d2 = d1;
   do {
      d2 >>= NTL_NBITS;
      i++;
   }
   while (d2 > 0);

   if (MustAlloc(a, i)) {
      zsetlength(&a, i);
      *aa = a;
   }

   i = 0;
   a[1] = 0;
   while (d1 > 0) {
      a[++i] = d1 & NTL_RADIXM;
      d1 >>= NTL_NBITS;
   }
   if (i > 0)
      a[0] = i;
   else
      a[0] = 1;

   if (anegative)
      a[0] = (-a[0]);
}

long ztoint(verylong a)
{
   long d;
   long sa;

   if (!a)
      return (0);

   if ((sa = *a) < 0)
      sa = -sa;

   d = *(a += sa);
   while (--sa) {
      d <<= NTL_NBITS;
      d += *(--a);
   }

   if ((*(--a)) < 0)
      return (-d);
   return (d);
}


long zcompare(verylong a, verylong b)
{
   long sa;
   long sb;

   if (!a) {
      if (!b)
         return (0);
      if (b[0] < 0)
         return (1);
      if (b[0] > 1)
         return (-1);
      if (b[1])
         return (-1);
      return (0);
   }
   if (!b) {
      if (a[0] < 0)
         return (-1);
      if (a[0] > 1)
         return (1);
      if (a[1])
         return (1);
      return (0);
   }

   if ((sa = *a) > (sb = *b))
      return (1);
   if (sa < sb)
      return (-1);
   if (sa < 0)
      sa = (-sa);
   a += sa;
   b += sa;
   for (; sa; sa--) {
      if (*a > *b) {
         if (sb < 0)
            return (-1);
         return (1);
      }
      if (*a-- < *b--) {
         if (sb < 0)
            return (1);
         return (-1);
      }
   }
   return (0);
}

void znegate(verylong *aa)
{
   verylong a = *aa;

   if (!a)
      return;
   if (a[1] || a[0] != 1)
      a[0] = (-a[0]);
}

void zsadd(verylong a, long d, verylong *b)
{
   verylong x = 0;

   zintoz(d, &x);
   zadd(a, x, b);
}


void
zadd(verylong a, verylong b, verylong *cc)
{
   long sa;
   long sb;
   long anegative;

   if (!a) {
      if (b)
         zcopy(b, cc);
      else
         zzero(cc);
      return;
   }

   if (!b) {
      zcopy(a, cc);
      return;
   }

   if ((anegative = ((sa = a[0]) < 0)) == ((sb = b[0]) < 0)) {
      /* signs a and b are the same */
      verylong pc, c;
      long carry;
      long i;
      long maxab;

      if (anegative) {
         sa = -sa;
         sb = -sb;
      }

      if (sa < sb) {
         i = sa;
         maxab = sb;
      }
      else {
         i = sb;
         maxab = sa;
      }

      c = *cc;
      if (MustAlloc(c, maxab+1)) {
         verylong c1 = c;
         zsetlength(&c1, maxab + 1);
         if (a == c) a = c1;
         if (b == c) b = c1;
         *cc = c = c1;
      }

      pc = c;
      carry = 0;

      do {
         long t = (*(++a)) + (*(++b)) + carry;
         carry = t >> NTL_NBITS;
         *(++pc) = t & NTL_RADIXM;
         i--;
      } while (i);

      i = sa-sb;
      if (!i) goto i_exit;

      if (i < 0) {
         i = -i;
         a = b;
      }

      if (!carry) goto carry_exit;

      for (;;) {
         long t = (*(++a)) + 1;
         carry = t >> NTL_NBITS;
         *(++pc) = t & NTL_RADIXM;
         i--;
         if (!i) goto i_exit;
         if (!carry) goto carry_exit;
      }

      i_exit:
      if (carry) {
         *(++pc) = 1;
         maxab++;
      }
      *c = anegative ? -maxab : maxab;
      return;

      carry_exit:
      if (pc != a) {
         do {
            *(++pc) = *(++a);
            i--;
         } while (i);
      }
      *c = anegative ? -maxab : maxab;
   }
   else {
      /* signs a and b are different...use zsub */

      verylong c = *cc;

      if (anegative) {
         a[0] = -sa;
         zsub(b, a, cc);
         if (a != c) a[0] = sa;
      }
      else {
         b[0] = -sb;
         zsub(a, b, cc);
         if (b != c) b[0] = sb;
      }
   }
}

void
zsub(verylong a, verylong b, verylong *cc)
{
   long sa;
   long sb;
   long anegative;

   if (!b) {
      if (a)
         zcopy(a, cc);
      else
         zzero(cc);
      return;
   }

   if (!a) {
      zcopy(b, cc);
      znegate(cc);
      return;
   }

   if ((anegative = ((sa = a[0]) < 0)) == ((sb = b[0]) < 0)) {
      /* signs agree */

      long i, carry, *pc, *c;

      if (anegative) {
         sa = -sa;
         sb = -sb;
      }

      carry = sa - sb;
      if (!carry) {
         long *aa = a + sa;
         long *bb = b + sa;
         
         i = sa;
         while (i && !(carry = (*aa - *bb))) {
            aa--; bb--; i--;
         }
      }

      if (!carry) {
         zzero(cc);
         return;
      }

      if (carry < 0) {
         { long t = sa; sa = sb; sb = t; }
         { long *t = a; a = b; b = t; }
         anegative = !anegative;
      }

      c = *cc;
      if (MustAlloc(c, sa)) {
         verylong c1 = c;
         zsetlength(&c1, sa);
         if (b == c) b = c1;
         *cc = c = c1;
      }

      i = sb;
      carry = 0;
      pc = c;

      do {
#if (!NTL_ARITH_RIGHT_SHIFT)
         long t = (*(++a)) - (*(++b)) - carry;
         carry = (t < 0);
#else
         long t = (*(++a)) - (*(++b)) + carry;
         carry = t >> NTL_NBITS;
#endif
         *(++pc) = t & NTL_RADIXM;
         i--;
      } while (i);

      i = sa-sb;
      while (carry) {
         long t = (*(++a)) - 1;
#if (!NTL_ARITH_RIGHT_SHIFT)
         carry = (t < 0);
#else
         carry = t >> NTL_NBITS;
#endif
         *(++pc) = t & NTL_RADIXM;
         i--;
      }

      if (i) {
         if (pc != a) {
            do {
               *(++pc) = *(++a);
               i--;
            } while (i);
         }
      }
      else {
         while (sa > 1 && *pc == 0) { sa--; pc--; } 
      }

      if (anegative) sa = -sa;
      *c = sa;
   }
   else {
      /* signs of a and b are different...use zadd */

      verylong c = *cc;

      if (anegative) {
         a[0] = -sa;
         zadd(a, b, cc);
         if (a != c) a[0] = sa;
         c = *cc;
         c[0] = -c[0];
      }
      else {
         b[0] = -sb;
         zadd(a, b, cc);
         if (b != c) b[0] = sb;
      }
   }
}


void
zsmul(verylong a, long d, verylong *bb)
{
   long sa;
   long anegative, bnegative;
   verylong b;


   if (d == 2) {
      z2mul(a, bb);
      return;
   }


   if ((d >= NTL_RADIX) || (d <= -NTL_RADIX)) {
      verylong x = 0;
      zintoz(d,&x);
      zmul(a, x, bb);
      return;
   }

   if (!a) {
      zzero(bb);
      return;
   }

   if (!d) {
      zzero(bb);
      return;
   }


   anegative = 0;
   bnegative = 0;

   if ((sa = a[0]) < 0) {
      anegative = 1;
      a[0] = sa = (-sa);
      if (d < 0)
         d = (-d);
      else
         bnegative = 1;
   }
   else if (bnegative = (d < 0))
      d = (-d);

   b = *bb;

   if (MustAlloc(b, sa + 1)) {
      zsetlength(&b, sa + 1);
      if (a == *bb) a = b;
      *bb = b;
   }

   zsxmul(d, b+1, a);

   sa++;
   while ((sa > 1) && (!(b[sa])))
      sa--;
   b[0] = sa;

   if (bnegative && (b[1] || b[0] != 1))
      b[0] = (-b[0]);

   if (anegative && a != b)
      a[0] = -a[0];
}

void zsubpos(verylong a, verylong b, verylong *cc)
{
   long sa;
   long sb;

   long *c, *pc;
   long i, carry;

   if (!b) {
      if (a)
         zcopy(a, cc);
      else
         zzero(cc);
      return;
   }

   if (!a) {
      zzero(cc);
      return;
   }

   sa = a[0];
   sb = b[0];

   c = *cc;
   if (MustAlloc(c, sa)) {
      verylong c1 = c;
      zsetlength(&c1, sa);
      if (b == c) b = c1;
      *cc = c = c1;
   }

   i = sb;
   carry = 0;
   pc = c;

   while (i) {
#if (!NTL_ARITH_RIGHT_SHIFT)
      long t = (*(++a)) - (*(++b)) - carry;
      carry = (t < 0);
#else
      long t = (*(++a)) - (*(++b)) + carry;
      carry = t >> NTL_NBITS;
#endif
      *(++pc) = t & NTL_RADIXM;
      i--;
   }

   i = sa-sb;
   while (carry) {
      long t = (*(++a)) - 1;
#if (!NTL_ARITH_RIGHT_SHIFT)
      carry = (t < 0);
#else
      carry = t >> NTL_NBITS;
#endif
      *(++pc) = t & NTL_RADIXM;
      i--;
   }

   if (i) {
      if (pc != a) {
         do {
            *(++pc) = *(++a);
            i--;
         } while (i);
      }
   }
   else {
      while (sa > 1 && *pc == 0) { sa--; pc--; } 
   }

   *c = sa;
}



static long *kmem = 0;     /* globals for Karatsuba */
static long max_kmem = 0;


/* These cross-over points were estimated using
   a Sparc-10, a Sparc-20, and a Pentium-90. */

#if (defined(NTL_SINGLE_MUL))
#define KARX (18)
#else
#define KARX (16) 
#endif

/* Auxilliary routines for Karatsuba */


static
void kar_fold(long *T, long *b, long hsa)
{
   long sb, *p2, *p3, i, carry;

   sb = *b;
   p2 = b + hsa;
   p3 = T;
   carry = 0;

   for (i = sb-hsa; i>0; i--) {
      long t = (*(++b)) + (*(++p2)) + carry;
      carry = t >> NTL_NBITS;
      *(++p3) = t & NTL_RADIXM;
   }

   for (i = (hsa << 1) - sb; i>0; i--) {
      long t = (*(++b)) + carry;
      carry = t >> NTL_NBITS;
      *(++p3) = t & NTL_RADIXM;
   }

   if (carry) {
      *(++p3) = carry;
      *T = hsa + 1;
   }
   else
      *T = hsa;
}

static
void kar_sub(long *T, long *c)
{
   long i, carry;

   i = *c;
   carry = 0;

   while (i>0) {
#if (!NTL_ARITH_RIGHT_SHIFT)
      long t = (*(++T)) - (*(++c)) - carry;
      carry = (t < 0);
#else
      long t = (*(++T)) - (*(++c)) + carry;
      carry = t >> NTL_NBITS;
#endif
      *T = t & NTL_RADIXM;
      i--;
   }

   while (carry) {
      long t = (*(++T)) - 1;
#if (!NTL_ARITH_RIGHT_SHIFT)
      carry = (t < 0);
#else
      carry = t >> NTL_NBITS;
#endif
      *T = t & NTL_RADIXM;
   }
}


static
void kar_add(long *c, long *T, long hsa)
{
   long i, carry;

   c += hsa;
   i = *T;
   while (T[i] == 0 && i > 0) i--;
   carry = 0;

   while (i>0) {
      long t = (*(++c)) + (*(++T)) + carry;
      carry = t >> NTL_NBITS;
      *c = t & NTL_RADIXM;
      i--;
   }

   while (carry) {
      long t = (*(++c)) + 1;
      carry = t >> NTL_NBITS;
      *c = t & NTL_RADIXM;
   }
}

static
void kar_fix(long *c, long *T, long hsa)
{
   long i, carry, s;

   s = *T;

   i = hsa;
   while (i>0) {
      *(++c) = *(++T);
      i--;
   }


   i = s - hsa;
   carry = 0;

   while (i > 0) {
      long t = (*(++c)) + (*(++T)) + carry;
      carry = t >> NTL_NBITS;
      *c = t & NTL_RADIXM;
      i--;
   }

   while (carry) {
      long t = (*(++c)) + 1;
      carry = t >> NTL_NBITS;
      *c = t & NTL_RADIXM;
   }
}
      
   

static
void kar_mul(long *c, long *a, long *b, long *stk)
{
   long sa, sb, sc;

   if (*a < *b) { long *t = a; a = b; b = t; }

   sa = *a;
   sb = *b;
   sc = sa + sb;

   if (sb < KARX) {
      /* classic algorithm */

      long *pc, i, *pb;

      pc = c;
      for (i = sc; i; i--) {
         pc++;
         *pc = 0;
      }
   
      pc = c;
      pb = b;
      for (i = sb; i; i--) {
         pb++;
         pc++;
         zaddmul(*pb, pc, a);
      }
   }
   else {
      long hsa = (sa + 1) >> 1;

      if (hsa < sb) {
         /* normal case */

         long *T1, *T2, *T3;

         /* allocate space */

         T1 = stk;  stk += hsa + 2;
         T2 = stk;  stk += hsa + 2;
         T3 = stk;  stk += (hsa << 1) + 3;

         if (stk-kmem > max_kmem) zhalt("internal error: kmem overflow");

         /* compute T1 = a_lo + a_hi */

         kar_fold(T1, a, hsa);

         /* compute T2 = b_lo + b_hi */

         kar_fold(T2, b, hsa);
         
         /* recursively compute T3 = T1 * T2 */

         kar_mul(T3, T1, T2, stk);

         /* recursively compute a_hi * b_hi into high part of c */
         /* and subtract from T3 */

         {
            long olda, oldb;

            olda = a[hsa];  a[hsa] = sa-hsa;
            oldb = b[hsa];  b[hsa] = sb-hsa;

            kar_mul(c + (hsa << 1), a + hsa, b + hsa, stk);
            kar_sub(T3, c + (hsa << 1));

            a[hsa] = olda;
            b[hsa] = oldb;
         }

         /* recursively compute a_lo*b_lo into low part of c */
         /* and subtract from T3 */

         *a = hsa;
         *b = hsa;

         kar_mul(c, a, b, stk);
         kar_sub(T3, c);

         *a = sa;
         *b = sb;

         /* finally, add T3 * NTL_RADIX^{hsa} to c */

         kar_add(c, T3, hsa);
      }
      else {
         /* degenerate case */

         long *T;
         
         T = stk;  stk += sb + hsa + 1;

         if (stk-kmem > max_kmem) zhalt("internal error: kmem overflow");

         /* recursively compute b*a_hi into high part of c */
         {
            long olda;

            olda = a[hsa];  a[hsa] = sa-hsa;
            kar_mul(c + hsa, a + hsa, b, stk);
            a[hsa] = olda;
         }

         /* recursively compute b*a_lo into T */

         *a = hsa;
         kar_mul(T, a, b, stk);
         *a = sa;

         /* fix-up result */

         kar_fix(c, T, hsa);
      }
   }

   /* normalize result */

   while (c[sc] == 0 && sc > 1) sc--;
   *c = sc;
}



#if (defined(NTL_SINGLE_MUL))
#define KARSX (36)
#else
#define KARSX (32) 
#endif

static
void kar_sq(long *c, long *a, long *stk)
{
   long sa, sc;

   sa = *a;
   sc = sa << 1;

   if (sa < KARSX) {
      /* classic algorithm */

      long carry, i, j, *pc;

      pc = c;
      for (i = sc; i; i--) {
         pc++;
         *pc = 0;
      }

      carry = 0;
      i = 0;
      for (j = 1; j <= sa; j++) {
         unsigned long uncar;
         long t;

         i += 2;
         uncar = ((unsigned long) carry) + (((unsigned long) c[i-1]) << 1);
         t = uncar & NTL_RADIXM;
         zaddmulpsq(t, a[j], carry);
         c[i-1] = t;
         zaddmulsq(sa-j, c+i, a+j);
         uncar = (c[i] << 1) + (uncar >> NTL_NBITS);
         uncar += carry;
         carry = uncar >> NTL_NBITS;
         c[i] = uncar & NTL_RADIXM;
      }

   }
   else {
      long hsa = (sa + 1) >> 1;
      long *T1, *T2, olda;

      T1 = stk;  stk += hsa + 2;
      T2 = stk;  stk += (hsa << 1) + 3;

      if (stk-kmem > max_kmem) zhalt("internal error: kmem overflow");

      kar_fold(T1, a, hsa);
      kar_sq(T2, T1, stk);

      olda = a[hsa];  a[hsa] = sa - hsa;
      kar_sq(c + (hsa << 1), a + hsa, stk);
      kar_sub(T2, c + (hsa << 1));
      a[hsa] = olda;

      *a = hsa;
      kar_sq(c, a, stk);
      kar_sub(T2, c);
      *a = sa;

      kar_add(c, T2, hsa);
   }

   while (c[sc] == 0 && sc > 1) sc--;
   *c = sc;
}

void zmul(verylong a, verylong b, verylong *cc)
{
   verylong mem = 0;
   verylong c = *cc;

   if (!a || (a[0] == 1 && a[1] == 0) || !b || (b[0] == 1 && b[1] == 0)) {
      zzero(cc);
      return;
   }

   if (a == b) {
      if (a == c) {
         zcopy(a, &mem);
         a = mem;
      }

      zsq(a, cc);
   }
   else {
      long aneg, bneg, sc;

      if (a == c) {
         zcopy(a, &mem);
         a = mem;
      }
      else if (b == c) {
         zcopy(b, &mem);
         b = mem;
      }

      if (*a < 0) {
         *a = - *a;
         aneg = 1;
      }
      else
         aneg = 0;

      if (*b < 0) {
         *b = - *b;
         bneg = 1;
      }
      else
         bneg = 0;

      sc = *a + *b;
      if (MustAlloc(c, sc)) {
         zsetlength(&c, sc);
         *cc = c;
      }

      if (*a < KARX || *b < KARX) {
         /* classic algorithm */

         long i, *pc;

         pc = c;
         for (i = sc; i; i--) {
            pc++;
            *pc = 0;
         }
   
         pc = c;
   
         if (*a >= *b) {
            long *pb = b;
            for (i = *pb; i; i--) {
               pb++;
               pc++;
               zaddmul(*pb, pc, a);
            }
         }
         else {
            long *pa = a;
            for (i = *pa; i; i--) {
               pa++; 
               pc++;
               zaddmul(*pa, pc, b);
            }
         }
   
         while (c[sc] == 0 && sc > 1) sc--;
         if (aneg != bneg && (sc != 1 || c[1] != 0)) sc = -sc;
         c[0] = sc;
         if (aneg) *a = - *a;
         if (bneg) *b = - *b;
      }
      else {
         /* karatsuba */

         long n, hn, sp;

         if (*a < *b)
            n = *b;
         else
            n = *a;

         sp = 0;
         do {
            hn = (n + 1) >> 1;
            sp += (hn << 2) + 7;
            n = hn+1;
         } while (n >= KARX);

         if (sp > max_kmem) {
            if (max_kmem == 0) 
               kmem = (long *) malloc(sp * sizeof(long));
            else
               kmem = (long *) realloc(kmem, sp * sizeof(long));

            max_kmem = sp;
            if (!kmem) zhalt("out of memory in karatsuba");
         }

         kar_mul(c, a, b, kmem);
         if (aneg != bneg && (*c != 1 || c[1] != 0)) *c = - *c;
         if (aneg) *a = - *a;
         if (bneg) *b = - *b;
      }
   }
}

void zsq(verylong a, verylong *cc)
{
   verylong mem = 0;
   verylong c = *cc;
   long aneg, sc;

   if (!a) {
      zzero(cc);
      return;
   }

   if (a == c) {
      zcopy(a, &mem);
      a = mem;
   }

   if (*a < 0) {
      *a = - *a;
      aneg = 1;
   }
   else
      aneg = 0;

   sc = (*a) << 1;
   if (MustAlloc(c, sc)) {
      zsetlength(&c, sc);
      *cc = c;
   }

   if (*a < KARSX) {
      /* classic algorithm */

      long carry, i, j, *pc, sa;

      pc = c;
      for (i = sc; i; i--) {
         pc++;
         *pc = 0;
      }

      carry = 0;
      i = 0;
      sa = *a;
      for (j = 1; j <= sa; j++) {
         unsigned long uncar;
         long t;

         i += 2;
         uncar = ((unsigned long) carry) + (((unsigned long) c[i-1]) << 1);
         t = uncar & NTL_RADIXM;
         zaddmulpsq(t, a[j], carry);
         c[i-1] = t;
         zaddmulsq(sa-j, c+i, a+j);
         uncar = (c[i] << 1) + (uncar >> NTL_NBITS);
         uncar += carry;
         carry = uncar >> NTL_NBITS;
         c[i] = uncar & NTL_RADIXM;
      }


      while (c[sc] == 0 && sc > 1) sc--;
      c[0] = sc;
      if (aneg) *a = - *a;
   }
   else {
      /* karatsuba */

      long n, hn, sp;

      n = *a;

      sp = 0;
      do {
         hn = (n + 1) >> 1;
         sp += hn + hn + hn + 5;
         n = hn+1;
      } while (n >= KARSX);

      if (sp > max_kmem) {
         if (max_kmem == 0) 
            kmem = (long *) malloc(sp * sizeof(long));
         else
            kmem = (long *) realloc(kmem, sp * sizeof(long));

         max_kmem = sp;
         if (!kmem) zhalt("out of memory in karatsuba");
      }

      kar_sq(c, a, kmem);
      if (aneg) *a = - *a;
   }
}

long zsdiv(verylong a, long d, verylong *bb)
{
   long sa;
   verylong b = *bb;

   if (!d) {
      zhalt("division by zero in zsdiv");
   }

   if (!a) {
      zzero(bb);
      return (0);
   }


   if (d == 2) {
      long is_odd = a[1] & 1;
      long fix = (a[0] < 0) & is_odd;
      zrshift(a, 1, bb);
      if (fix) zsadd(*bb, -1, bb);
      return is_odd;
   } 


   if ((sa = a[0]) < 0)
      sa = (-sa);
   zsetlength(&b, sa);
   if (a == *bb) a = b;
   *bb = b;
   if ((d >= NTL_RADIX) || (d <= -NTL_RADIX)) {
      verylong zd = 0;
      verylong zb = 0;

      zintoz(d, &zb);
      zdiv(a, zb, &b, &zd);
      *bb = b;
      return (ztoint(zd));
   }
   else {
      long den = d;
      double deninv;
      long carry = 0;
      long i;
      long flag = (*a < 0 ? 2 : 0) | (den < 0 ? 1 : 0);

      if (den < 0)
         den = -den;
      deninv = (double)1/den;
      if (a[sa] < den && sa > 1)
         carry = a[sa--];
      for (i = sa; i; i--) {
         zdiv21(carry, a[i], den, deninv, b[i]);
      }
      while ((sa > 1) && (!(b[sa])))
         sa--;
      b[0] = sa;
      if (flag) {
         if (flag <= 2) {
            if (!carry)
               znegate(&b);
            else {
               zsadd(b, 1, &b);
               b[0] = -b[0];
               if (flag == 1)
                  carry = carry - den;
               else
                  carry = den - carry;
               *bb = b;
            }
         }
         else
            carry = -carry;
      }
      return (carry);
   }
}

long zsmod(verylong a, long d)
{
   long sa;

   if (!a) {
      return (0);
   }

   if (d == 2) return (a[1] & 1);

   if (!d) {
      zhalt("division by zero in zsdiv");
   }

   if ((sa = a[0]) < 0)
      sa = (-sa);
   if ((d >= NTL_RADIX) || (d <= -NTL_RADIX)) {
      verylong zd = 0;
      verylong zb = 0;

      zintoz(d, &zb);
      zmod(a, zb, &zd);
      return (ztoint(zd));
   }
   else {
      long den = d;
      double deninv;
      long carry = 0;
      long i;
      long flag = (*a < 0 ? 2 : 0) | (den < 0 ? 1 : 0);

      if (den < 0)
         den = -den;
      deninv = (double)1/den;
      if (a[sa] < den && sa > 1)
         carry = a[sa--];
      for (i = sa; i; i--) {
         zrem21(carry, a[i], den, deninv);
      }
      if (flag) {
         if (flag <= 2) {
            if (carry) {
               if (flag == 1)
                  carry = carry - den;
               else
                  carry = den - carry;
            }
         }
         else
            carry = -carry;
      }
      return (carry);
   }
}

void zmultirem(verylong a, long n, long* dd, long *rr)
{
   long j;
   long sa;

   if (!a || (a[0] == 1 && a[1] == 0)) {
      for (j = 0; j < n; j++) rr[j] = 0;
      return;
   }

   sa = a[0];

   for (j = 0; j < n; j++) {
      long den = dd[j];
      double deninv;
      long carry = 0;
      long i;
      long lsa = sa;

      deninv = (double)1/den;
      if (a[lsa] < den && lsa > 1)
         carry = a[lsa--];
      for (i = lsa; i; i--) {
         zrem21(carry, a[i], den, deninv);
      }
      rr[j] = carry;
   }
}


#if (defined(NTL_SINGLE_MUL))

#define REDUCE_CROSS 8

void zmultirem2(verylong a, long n, long* dd, double **ttbl, long *rr)

{
   long d; 
   double *tbl;
   long ac0, ac1, ac2;
   long *ap; 
   double *tp;
   long sa, i, j, k;
   unsigned long pc0, pc1, lo_wd, hi_wd;
   long carry;
   double xx, yy, dinv;

   if (!a ||  a[0] < REDUCE_CROSS || a[0] >= NTL_RADIX) {
      zmultirem(a, n, dd, rr);
      return;
   }

   sa = a[0];

   for (i = 0; i < n; i++) {
      d = dd[i];
      tbl = ttbl[i];

      ac0 = a[1];
      ac1 = 0;
      ac2 = 0;

      ap = &a[2];
      tp = &tbl[1];

      k = sa-1;
      while (k) {
         j = k;
         if (j > 64) j = 64;
         k -= j;

         pc0 = 0;
         pc1 = 0;

         xx =  ((double) (*ap++))* (*tp++) + DENORMALIZE;
         for (; j > 1; j--) {
            yy =  ((double) (*ap++))*(*tp++) 
                  + DENORMALIZE;
            NTL_FetchHiLo(hi_wd, lo_wd, xx);
            pc0 +=  lo_wd & 0x3FFFFFF;
            pc1 += ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
            xx = yy;
         }
         NTL_FetchHiLo(hi_wd, lo_wd, xx);
         pc0 +=  lo_wd & 0x3FFFFFF;
         pc1 += ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
         pc1 += (pc0 >> 26);

         ac0 += (pc0 & 0x3FFFFFF);
         if (ac0 > NTL_RADIX) {
            ac0 -= NTL_RADIX;
            carry = 1;
         }
         else
            carry = 0;

         ac1 += carry + (pc1 & 0x3FFFFFF);
         if (ac1 > NTL_RADIX) {
            ac1 -= NTL_RADIX;
            carry = 1;
         }
         else
            carry = 0;

         ac2 += carry + (pc1 >> 26);
      }

      carry = 0;
      dinv = ((double) 1)/((double) d);
      if (ac2 >= d) {
         zrem21(carry, ac2, d, dinv);
      }
      else
         carry = ac2;

      zrem21(carry, ac1, d, dinv);
      zrem21(carry, ac0, d, dinv);

      rr[i] = carry;
   }
}


#elif (defined(NTL_TBL_REM))

void zmultirem3(verylong a, long n, long* dd, long **ttbl, long *rr)
{
   long sa, i, d, *tbl, ac0, ac1, ac2, *ap, *tp, k, t, carry;
   double dinv;

   if (!a || a[0] < 8 || a[0] >= NTL_RADIX) {
      zmultirem(a, n, dd, rr);
      return;
   }

   sa = a[0];
 
   for (i = 0; i < n; i++) {
      d = dd[i];
      tbl = ttbl[i];
      ac0 = a[1];
      ac1 = 0;
      ac2 = 0;
      ap = &a[2];
      tp = &tbl[1];
      k = sa-1;
      
      while (k) {
         zxmulp(t, *ap, *tp, ac0);
         ac1 += ac0;
         ac2 += ac1 >> NTL_NBITS;
         ac1 &= NTL_RADIXM;
         ac0 = t;
         k--;
         ap++;
         tp++;
      }

      carry = 0;
      dinv = ((double) 1)/((double) d);
      if (ac2 >= d) {
         zrem21(carry, ac2, d, dinv);
      }
      else
         carry = ac2;

      zrem21(carry, ac1, d, dinv);
      zrem21(carry, ac0, d, dinv);

      rr[i] = carry;
   }
}

#endif




long zsfastrem(verylong a, long d)
/* assumes a >= 0, and 0 < d < NTL_RADIX, i=0..n-1 */
/* computes a % d */

{
   long sa;

   if (!a || (a[0] == 1 && a[1] == 0)) {
      return 0;
   }

   sa = a[0];

   {
      long den = d;
      double deninv = ((double)1)/((double)den);
      long carry = 0;
      long i;
      long lsa = sa;

      if (a[lsa] < den && lsa > 1)
         carry = a[lsa--];
      for (i = lsa; i; i--) {
         zrem21(carry, a[i], den, deninv);
      }
      return carry;
   }
}



void zdiv(verylong a, verylong b, verylong *qq, verylong *rr)
{
   long sa, sb, sq, i;
   long sign;
   long q1;
   long *rp;
   double btopinv, aux;
   verylong q=0, r=0;

   if (!b || (((sb=b[0]) == 1) && (!b[1]))) {
      zhalt("division by zero in zdiv");
   }

   if (!a || (((sa=a[0]) == 1) && (!a[1]))) {
      zzero(qq);
      if (rr) zzero(rr);
      return;
   }


   if (sb == 1) {
      long t1 = zsdiv(a, b[1], qq);
      if (rr) zintoz(t1, rr);
      return;
   }

   if (sb == -1) {
      long t1 = zsdiv(a, -b[1], qq);
      if (rr) zintoz(t1, rr);
      return;
   }

   sign = 0;
   if (sa < 0) {
      a[0] = sa = -sa;
      sign = 2;
   }

   if (sb < 0) {
      b[0] = sb = -sb;
      sign |= 1;
   }


   sq = sa-sb+1;

   if (sq <= 0) {
      zcopy(a, &r);
      zzero(&q);
      goto done;
   }

   zsetlength(&q, sq);
   zsetlength(&r, sa+1);
 
   zcopy(a, &r);
   rp = &r[sa+1];
   *rp = 0;

   r[0] = 0; /* this streamlines the last evaluation of aux */

   btopinv = b[sb]*NTL_FRADIX + b[sb-1];
   if (sb > 2) 
      btopinv = NTL_FRADIX / (btopinv*NTL_FRADIX + b[sb-2]);
   else
      btopinv = 1.0 / btopinv;

   
   aux = btopinv*(rp[-1]*NTL_FRADIX + rp[-2]);
   if (aux >= NTL_FRADIX)
      aux = NTL_FRADIX-1;

   for (i = sq; i >= 1; i--, rp--) {
      q1 = (long) aux;
      if (q1) {
         zsubmul(q1, &r[i], b);
      }

      while (rp[0] < 0) {
         zaddmulone(&r[i], b);
         q1--;
      }

      while (rp[0] > 0) {
         zsubmulone(&r[i], b);
         q1++;
      }

      aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
      while (aux > NTL_FRADIX - 8) {
         /* q1 might be too small */
         if (aux >= NTL_FRADIX)
            aux = NTL_FRADIX-1;


         zsubmulone(&r[i], b);
         if (rp[0] < 0) {
            /* oops...false alarm! */
            zaddmulone(&r[i], b);
            break;
         }
         else {
            q1++;
            aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
         }
      }

      q[i] = q1;
   }

   while (sq > 1 && q[sq] == 0) sq--;
   q[0] = sq;

   i = sb;
   while (i > 1 && r[i] == 0) i--;
   r[0] = i;

done:
   if (sign)
   {
      if (sign <= 2)
      {
         if (!(r[1]) && (r[0] == 1))
            znegate(&q);
         else
         {
            zsadd(q, 1, &q);
            znegate(&q);
            if (sign == 1)
               zsub(r, b, &r);
            else
               zsub(b, r, &r);
         }
      }
      else
         znegate(&r);

      if (sign & 2)
         a[0] = -sa;

      if (sign & 1)
         b[0] = -sb;
   }

   zcopy(q, qq);
   if (rr) zcopy(r, rr);
}

void
zmod(verylong a, verylong b, verylong *rr)
{
   long sa, sb, sq, i;
   long sign;
   long q1;
   long *rp;
   double btopinv, aux;
   verylong r=0;

   if (!b || (((sb=b[0]) == 1) && (!b[1]))) {
      zhalt("division by zero in zdiv");
   }

   if (!a || (((sa=a[0]) == 1) && (!a[1]))) {
      zzero(rr);
      return;
   }


   if (sb == 1) {
      zintoz(zsmod(a, b[1]), rr);
      return;
   }

   if (sb == -1) {
      zintoz(zsmod(a, -b[1]), rr);
      return;
   }

   sign = 0;
   if (sa < 0) {
      a[0] = sa = -sa;
      sign = 2;
   }

   if (sb < 0) {
      b[0] = sb = -sb;
      sign |= 1;
   }


   sq = sa-sb+1;

   if (sq <= 0) {
      zcopy(a, &r);
      goto done;
   }

   zsetlength(&r, sa+1);
 
   zcopy(a, &r);
   rp = &r[sa+1];
   *rp = 0;

   r[0] = 0; /* this streamlines the last evaluation of aux */

   btopinv = b[sb]*NTL_FRADIX + b[sb-1];
   if (sb > 2) 
      btopinv = NTL_FRADIX / (btopinv*NTL_FRADIX + b[sb-2]);
   else
      btopinv = 1.0 / btopinv;

   
   aux = btopinv*(rp[-1]*NTL_FRADIX + rp[-2]);
   if (aux >= NTL_FRADIX)
      aux = NTL_FRADIX-1;

   for (i = sq; i >= 1; i--, rp--) {
      q1 = (long) aux;
      if (q1) {
         zsubmul(q1, &r[i], b);
      }

      while (rp[0] < 0) {
         zaddmulone(&r[i], b);
      }

      while (rp[0] > 0) {
         zsubmulone(&r[i], b);
      }

      aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
      while (aux > NTL_FRADIX - 8) {
         /* q1 might be too small */
         if (aux >= NTL_FRADIX)
            aux = NTL_FRADIX-1;


         zsubmulone(&r[i], b);
         if (rp[0] < 0) {
            /* oops...false alarm! */
            zaddmulone(&r[i], b);
            break;
         }
         else {
            aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
         }
      }
   }

   i = sb;
   while (i > 1 && r[i] == 0) i--;
   r[0] = i;

done:
   if (sign)
   {
      if (sign <= 2)
      {
         if (!(r[1]) && (r[0] == 1))
            /* no op */;
         else
         {
            if (sign == 1)
               zsub(r, b, &r);
            else
               zsub(b, r, &r);
         }
      }
      else
         znegate(&r);

      if (sign & 2)
         a[0] = -sa;

      if (sign & 1)
         b[0] = -sb;
   }

   zcopy(r, rr);
}

void
zquickmod(verylong *rr, verylong b)
{
   long sa, sb, sq, i;
   long q1;
   long *rp;
   double btopinv, aux;
   verylong r;

   sb = b[0];

   r = *rr;

   if (!r || (((sa=r[0]) == 1) && (!r[1]))) {
      zzero(rr);
      return;
   }


   if (sb == 1) {
      zintoz(zsmod(r, b[1]), rr);
      return;
   }

   sq = sa-sb+1;

   if (sq <= 0) {
      return;
   }

   zsetlength(rr, sa+1);
   r = *rr;
 
   rp = &r[sa+1];
   *rp = 0;

   r[0] = 0; /* this streamlines the last evaluation of aux */

   btopinv = b[sb]*NTL_FRADIX + b[sb-1];
   if (sb > 2) 
      btopinv = NTL_FRADIX / (btopinv*NTL_FRADIX + b[sb-2]);
   else
      btopinv = 1.0 / btopinv;

   
   aux = btopinv*(rp[-1]*NTL_FRADIX + rp[-2]);
   if (aux >= NTL_FRADIX)
      aux = NTL_FRADIX-1;

   for (i = sq; i >= 1; i--, rp--) {
      q1 = (long) aux;
      if (q1) {
         zsubmul(q1, &r[i], b);
      }

      while (rp[0] < 0) {
         zaddmulone(&r[i], b);
      }

      while (rp[0] > 0) {
         zsubmulone(&r[i], b);
      }

      aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
      while (aux > NTL_FRADIX - 8) {
         /* q1 might be too small */
         if (aux >= NTL_FRADIX)
            aux = NTL_FRADIX-1;


         zsubmulone(&r[i], b);
         if (rp[0] < 0) {
            /* oops...false alarm! */
            zaddmulone(&r[i], b);
            break;
         }
         else {
            aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
         }
      }
   }

   i = sb;
   while (i > 1 && r[i] == 0) i--;
   r[0] = i;
}


void
zaddmod(
	verylong a,
	verylong b,
	verylong n,
	verylong *c
	)
{
	zadd(a, b, c);
	if (zcompare(*c, n) >= 0)
		zsubpos(*c, n, c);
}

void
zsubmod(
	verylong a,
	verylong b,
	verylong n,
	verylong *c
	)
{
	verylong mem = 0;
	long cmp;

	if ((cmp=zcompare(a, b)) < 0) {
		zadd(n, a, &mem);
		zsubpos(mem, b, c);
	} else if (!cmp) 
		zzero(c);
	else 
		zsubpos(a, b, c);
}

void
zsmulmod(
	verylong a,
	long d,
	verylong n,
	verylong *c
	)
{
	verylong mem = 0;

	zsmul(a, d, &mem);
	zquickmod(&mem, n);
	zcopy(mem, c);
}

void
zmulmod(
	verylong a,
	verylong b,
	verylong n,
	verylong *c
	)
{
	verylong mem = 0;

	zmul(a, b, &mem);
	zquickmod(&mem, n);
	zcopy(mem, c);
}

void
zsqmod(
	verylong a,
	verylong n,
	verylong *c
	)
{
	verylong mem = 0;

	zsq(a, &mem);
	zquickmod(&mem, n);
	zcopy(mem, c);
}



void
zinvmod(
	verylong a,
	verylong n,
	verylong *c
	)
{
	if (zinv(a, n, c))
		zhalt("undefined inverse in zinvmod");
}


long 
zxxeucl(
   verylong ain,
   verylong nin,
   verylong *invv,
   verylong *uu
   )
{
   verylong a = 0;
   verylong n = 0;
   verylong q = 0;
   verylong w = 0;
   verylong x = 0;
   verylong y = 0;
   verylong z = 0;
   verylong inv = *invv;
   verylong u = *uu;
   long diff;
   long ilo;
   long sa;
   long sn;
   long temp;
   long e;
   long fast;
   long parity;
   long gotthem;
   verylong pin;
   verylong p;
   long i;
   long try11;
   long try12;
   long try21;
   long try22;
   long got11;
   long got12;
   long got21;
   long got22;
   double hi;
   double lo;
   double dt;
   double fhi, fhi1;
   double flo, flo1;
   double num;
   double den;
   double dirt;


   zsetlength(&a, (e = (ain[0] > nin[0] ? ain[0] : nin[0])));
   zsetlength(&n, e);
   zsetlength(&q, e);
   zsetlength(&w, e);
   zsetlength(&x, e);
   zsetlength(&y, e);
   zsetlength(&z, e);
   zsetlength(&inv, e);
   *invv = inv;
   zsetlength(&u, e);
   *uu = u;

   fhi1 = 1.0 + ((double) 32.0)/NTL_FDOUBLE_PRECISION;
   flo1 = 1.0 - ((double) 32.0)/NTL_FDOUBLE_PRECISION;

   fhi = 1.0 + ((double) 8.0)/NTL_FDOUBLE_PRECISION;
   flo = 1.0 - ((double) 8.0)/NTL_FDOUBLE_PRECISION;

   pin = &ain[0];
   p = &a[0];
   for (i = (*pin); i >= 0; i--)
      *p++ = *pin++;
   pin = &nin[0];
   p = &n[0];
   for (i = (*pin); i >= 0; i--)
      *p++ = *pin++;
   inv[0] = 1;
   inv[1] = 1;
   w[0] = 1;
   w[1] = 0;
   while (n[0] > 1 || n[1] > 0)
   {
      gotthem = 0;
      sa = a[0];
      sn = n[0];
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = a[0];
         p = &a[sa];
         num = (double) (*p) * NTL_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_FRADIX;
         if (sa > 2)
            num += (*(p - 1));
         sn = n[0];
         p = &n[sn];
         den = (double) (*p) * NTL_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_FRADIX;
         if (sn > 2)
            den += (*(p - 1));
         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_FRADIX;
            lo *= NTL_FRADIX;
         }
         try11 = 1;
         try12 = 0;
         try21 = 0;
         try22 = 1;
         parity = 1;
         fast = 1; 
         while (fast > 0)
         {
            parity = 1 - parity;
            if (hi >= NTL_FRADIX)
               fast = 0;
            else
            {
               ilo = (long)lo;
               dirt = hi - ilo;
               if (dirt <= 0 || !ilo || ilo < (long)hi)
                  fast = 0;
               else
               {
                  dt = lo-ilo;
                  lo = flo / dirt;
                  if (dt > 0)
                     hi = fhi / dt;
                  else
                     hi = NTL_FRADIX;
                  temp = try11;
                  try11 = try21;
                  if ((NTL_RADIX - temp) / ilo < try21)
                     fast = 0;
                  else
                     try21 = temp + ilo * try21;
                  temp = try12;
                  try12 = try22;
                  if ((NTL_RADIX - temp) / ilo < try22)
                     fast = 0;
                  else
                     try22 = temp + ilo * try22;
                  if ((fast > 0) && (parity > 0))
                  {
                     gotthem = 1;
                     got11 = try11;
                     got12 = try12;
                     got21 = try21;
                     got22 = try22;
                  }
               }
            }
         }
      }
      if (gotthem)
      {
         zsmul(inv, got11, &x);
         zsmul(w, got12, &y);
         zsmul(inv, got21, &z);
         zsmul(w, got22, &w);
         zadd(x, y, &inv);
         zadd(z, w, &w);
         zsmul(a, got11, &x);
         zsmul(n, got12, &y);
         zsmul(a, got21, &z);
         zsmul(n, got22, &n);
         zsub(x, y, &a);
         zsub(n, z, &n);
      }
      else
      {
         zdiv(a, n, &q, &a);
         zmul(q, w, &x);
         zadd(inv, x, &inv);
         if (a[0] > 1 || a[1] > 0)
         {
            zdiv(n, a, &q, &n);
            zmul(q, inv, &x);
            zadd(w, x, &w);
         }
         else
         {
            p = &a[0];
            pin = &n[0];
            for (i = (*pin); i >= 0; i--)
               *p++ = *pin++;
            n[0] = 1;
            n[1] = 0;
            zcopy(w, &inv);
            znegate(&inv);
         }
      }
   }

   if ((a[0] == 1) && (a[1] == 1))
      e = 0;
   else 
      e = 1;

   p = &u[0];
   pin = &a[0];
   for (i = (*pin); i >= 0; i--)
      *p++ = *pin++;
   *invv = inv;
   *uu = u;

   return (e);
}

long 
zinv(
	verylong ain,
	verylong nin,
	verylong *invv
	)
{
	verylong u = 0;
	verylong v = 0;
	long sgn;


	if (zscompare(nin, 1) <= 0) {
		zhalt("InvMod: second input <= 1");
	}

	sgn = zsign(ain);
	if (sgn < 0) {
		zhalt("InvMod: first input negative");
	}

	if (zcompare(ain, nin) >= 0) {
		zhalt("InvMod: first input too big");
	}


	if (sgn == 0) {
		zcopy(nin, invv);
		return 1;
	}

	if (!(zxxeucl(ain, nin, &v, &u))) {
		if (zsign(v) < 0) zadd(v, nin, &v);
		zcopy(v, invv);
		return 0;
	}

	zcopy(u, invv);
	return 1;
}

void
zexteucl(
	verylong aa,
	verylong *xa,
	verylong bb,
	verylong *xb,
	verylong *d
	)
{
	verylong modcon = 0;
	verylong a=0, b=0;
	long anegative = 0;
	long bnegative = 0;

	zcopy(aa, &a);
	zcopy(bb, &b);

	if (anegative = (a[0] < 0))
		a[0] = -a[0];
	if (bnegative = (b[0] < 0))
		b[0] = -b[0];

	if (!b[1] && (b[0] == 1))
	{
		zone(xa);
		zzero(xb);
		zcopy(a, d);
		goto done;
	}

	if (!a[1] && (a[0] == 1))
	{
		zzero(xa);
		zone(xb);
		zcopy(b, d);
		goto done;
	}

	zxxeucl(a, b, xa, d);
	zmul(a, *xa, xb);
	zsub(*d, *xb, xb);
	zdiv(*xb, b, xb, &modcon);

	if ((modcon[1]) || (modcon[0] != 1))
	{
		zhalt("non-zero remainder in zexteucl   BUG");
	}
done:
	if (anegative)
	{
		znegate(xa);
	}
	if (bnegative)
	{
		znegate(xb);
	}
}



/* I've adapted LIP's extended euclidean algorithm to
 * do rational reconstruction.  -- VJS.
 */


long 
zxxratrecon(
   verylong ain,
   verylong nin,
   verylong num_bound,
   verylong den_bound,
   verylong *num_out,
   verylong *den_out
   )
{
   verylong a = 0;
   verylong n = 0;
   verylong q = 0;
   verylong w = 0;
   verylong x = 0;
   verylong y = 0;
   verylong z = 0;
   verylong inv = 0;
   verylong u = 0;
   verylong a_bak = 0;
   verylong n_bak = 0;
   verylong inv_bak = 0;
   verylong w_bak = 0;

   verylong p;

   long diff;
   long ilo;
   long sa;
   long sn;
   long snum;
   long sden;
   long e;
   long fast;
   long temp;
   long parity;
   long gotthem;
   long try11;
   long try12;
   long try21;
   long try22;
   long got11;
   long got12;
   long got21;
   long got22;

   double hi;
   double lo;
   double dt;
   double fhi, fhi1;
   double flo, flo1;
   double num;
   double den;
   double dirt;

   if (zsign(num_bound) < 0)
      zhalt("rational reconstruction: bad numerator bound");

   if (!num_bound)
      snum = 1;
   else
      snum = num_bound[0];

   if (zsign(den_bound) <= 0)
      zhalt("rational reconstruction: bad denominator bound");

   sden = den_bound[0];

   if (zsign(nin) <= 0)
      zhalt("rational reconstruction: bad modulus");

   if (zsign(ain) < 0 || zcompare(ain, nin) >= 0)
      zhalt("rational reconstruction: bad residue");

      
   zsetlength(&a, (e = (ain[0] > nin[0] ? ain[0] : nin[0])));
   zsetlength(&n, e);
   zsetlength(&q, e);
   zsetlength(&w, e);
   zsetlength(&x, e);
   zsetlength(&y, e);
   zsetlength(&z, e);
   zsetlength(&inv, e);
   zsetlength(&u, e);
   zsetlength(&a_bak, e);
   zsetlength(&n_bak, e);
   zsetlength(&inv_bak, e);
   zsetlength(&w_bak, e);

   fhi1 = 1.0 + ((double) 32.0)/NTL_FDOUBLE_PRECISION;
   flo1 = 1.0 - ((double) 32.0)/NTL_FDOUBLE_PRECISION;

   fhi = 1.0 + ((double) 8.0)/NTL_FDOUBLE_PRECISION;
   flo = 1.0 - ((double) 8.0)/NTL_FDOUBLE_PRECISION;

   zcopy(ain, &a);
   zcopy(nin, &n);

   zone(&inv);
   zzero(&w);

   while (1)
   {
      if (w[0] >= sden && zcompare(w, den_bound) > 0) break;
      if (n[0] <= snum && zcompare(n, num_bound) <= 0) break;

      zcopy(a, &a_bak);
      zcopy(n, &n_bak);
      zcopy(w, &w_bak);
      zcopy(inv, &inv_bak);

      gotthem = 0;
      sa = a[0];
      sn = n[0];
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = a[0];
         p = &a[sa];
         num = (double) (*p) * NTL_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_FRADIX;
         if (sa > 2)
            num += (*(p - 1));
         sn = n[0];
         p = &n[sn];
         den = (double) (*p) * NTL_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_FRADIX;
         if (sn > 2)
            den += (*(p - 1));
         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_FRADIX;
            lo *= NTL_FRADIX;
         }
         try11 = 1;
         try12 = 0;
         try21 = 0;
         try22 = 1;
         parity = 1;
         fast = 1; 
         while (fast > 0)
         {
            parity = 1 - parity;
            if (hi >= NTL_FRADIX)
               fast = 0;
            else
            {
               ilo = (long)lo;
               dirt = hi - ilo;
               if (dirt <= 0 || !ilo || ilo < (long)hi)
                  fast = 0;
               else
               {
                  dt = lo-ilo;
                  lo = flo / dirt;
                  if (dt > 0)
                     hi = fhi / dt;
                  else
                     hi = NTL_FRADIX;
                  temp = try11;
                  try11 = try21;
                  if ((NTL_RADIX - temp) / ilo < try21)
                     fast = 0;
                  else
                     try21 = temp + ilo * try21;
                  temp = try12;
                  try12 = try22;
                  if ((NTL_RADIX - temp) / ilo < try22)
                     fast = 0;
                  else
                     try22 = temp + ilo * try22;
                  if ((fast > 0) && (parity > 0))
                  {
                     gotthem = 1;
                     got11 = try11;
                     got12 = try12;
                     got21 = try21;
                     got22 = try22;
                  }
               }
            }
         }
      }
      if (gotthem)
      {
         zsmul(inv, got11, &x);
         zsmul(w, got12, &y);
         zsmul(inv, got21, &z);
         zsmul(w, got22, &w);
         zadd(x, y, &inv);
         zadd(z, w, &w);
         zsmul(a, got11, &x);
         zsmul(n, got12, &y);
         zsmul(a, got21, &z);
         zsmul(n, got22, &n);
         zsub(x, y, &a);
         zsub(n, z, &n);
      }
      else
      {
         zdiv(a, n, &q, &a);
         zmul(q, w, &x);
         zadd(inv, x, &inv);
         if (a[0] > 1 || a[1] > 0)
         {
            zdiv(n, a, &q, &n);
            zmul(q, inv, &x);
            zadd(w, x, &w);
         }
         else
         {
            break;
         }
      }
   }

   zcopy(a_bak, &a);
   zcopy(n_bak, &n);
   zcopy(w_bak, &w);
   zcopy(inv_bak, &inv);

   znegate(&w);

   while (1)
   {
      sa = w[0];
      if (sa < 0) w[0] = -sa;
      if (w[0] >= sden && zcompare(w, den_bound) > 0) return 0;
      w[0] = sa;

      if (n[0] <= snum && zcompare(n, num_bound) <= 0) break;
      
      fast = 0;
      sa = a[0];
      sn = n[0];
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = a[0];
         p = &a[sa];
         num = (double) (*p) * NTL_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_FRADIX;
         if (sa > 2)
            num += (*(p - 1));
         sn = n[0];
         p = &n[sn];
         den = (double) (*p) * NTL_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_FRADIX;
         if (sn > 2)
            den += (*(p - 1));
         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_FRADIX;
            lo *= NTL_FRADIX;
         }
         if (hi < NTL_FRADIX)
         {
            ilo = (long)lo;
            if (ilo == (long)hi)
               fast = 1;
         }
      }

      if (fast) 
      {
         if (ilo != 0) {
            if (ilo == 1) {
               zsub(inv, w, &inv);
               zsubpos(a, n, &a);
            }
            else if (ilo == 2) {
               z2mul(w, &x);
               zsub(inv, x, &inv);
               z2mul(n, &x);
               zsubpos(a, x, &a);
            }
            else if (ilo ==3) {
               z2mul(w, &x);
               zadd(w, x, &x);
               zsub(inv, x, &inv);
               z2mul(n, &x);
               zadd(n, x, &x);
               zsubpos(a, x, &a);
            }
            else if (ilo == 4) {
               zlshift(w, 2, &x);
               zsub(inv, x, &inv);
               zlshift(n, 2, &x);
               zsubpos(a, x, &a);
            }
            else {
               zsmul(w, ilo, &x);
               zsub(inv, x, &inv);
               zsmul(n, ilo, &x);
               zsubpos(a, x, &a);
            }
         }
      }
      else {
         zdiv(a, n, &q, &a);
         zmul(q, w, &x);
         zsub(inv, x, &inv);
      }

      zswap(&a, &n);
      zswap(&inv, &w);
   }

   if (zsign(w) < 0) {
      znegate(&w);
      znegate(&n);
   }

   zcopy(n, num_out);
   zcopy(w, den_out);

   return 1;
}


void
zsexpmod(
	verylong a,
	long e,
	verylong n,
	verylong *bb
	)
{
	verylong le = 0;

	zintoz(e, &le);
	zexpmod(a, le, n, bb);
}

void
zexpmod(
	verylong a,
	verylong e,
	verylong n,
	verylong *bb
	)
{
	long i;
	long j;
	long k = 0;
	verylong res = 0;

	if (zsign(a) < 0 || zcompare(a, n) >= 0 || zscompare(n, 1) <= 0)
		zhalt("zexpmod: bad args");

	if (!e)
	{
		zone(bb);
		return;
	}

	if ((i = e[0]) < 0)
		zhalt("negative exponent in zexpmod");

	if ((i == 1) && (!(e[1]))) {
		zone(bb);
		return;
	}

	if (ziszero(a))
	{
		zzero(bb);
		return;
	}

	zcopy(a, &res);
	for (; i; i--)
	{
		j = e[i];
		if (!k)
		{
			k = 1;
			while ((k << 1) <= j)
				k <<= 1;
		}
		while (k >>= 1)
		{
			zsqmod(res, n, &res);
			if (j & k)
				zmulmod(a, res, n, &res);
		}
		k = NTL_RADIX;
	}
	zcopy(res, bb);
}

void
zexp(
	verylong a,
	long e,
	verylong *bb
	)
{
	long k;
	long len_a;
	verylong res = 0;

	if (!e)
	{
		zone(bb);
		return;
	}

	if (e < 0)
		zhalt("negative exponent in zexp");

	if (ziszero(a))
	{
		zzero(bb);
		return;
	}

	len_a = z2log(a);
	if (len_a > (NTL_MAX_LONG-(NTL_NBITS-1))/e)
		zhalt("overflow in zexp");

	zsetlength(&res, (len_a*e+NTL_NBITS-1)/NTL_NBITS);

	zcopy(a, &res);
	k = 1;
	while ((k << 1) <= e)
		k <<= 1;
	while (k >>= 1) {
		zsq(res, &res);
		if (e & k)
			zmul(a, res, &res);
	}

	zcopy(res, bb);
}

void
zexps(
	long a,
	long e,
	verylong *bb
	)
{
	long k;
	long len_a;
	verylong res = 0;

	if (!e)
	{
		zone(bb);
		return;
	}

	if (e < 0)
		zhalt("negative exponent in zexps");

	if (!a)
	{
		zzero(bb);
		return;
	}

	if (a >= NTL_RADIX || a <= -NTL_RADIX) {
		zintoz(a, &res);
		zexp(res, e, &res);
		return;
	}

	len_a = z2logs(a);
	if (len_a > (NTL_MAX_LONG-(NTL_NBITS-1))/e)
		zhalt("overflow in zexps");

	zsetlength(&res, (len_a*e+NTL_NBITS-1)/NTL_NBITS);

	zintoz(a, &res);
	k = 1;
	while ((k << 1) <= e)
		k <<= 1;
	while (k >>= 1) {
		zsq(res, &res);
		if (e & k)
			zsmul(res, a, &res);
	}

	zcopy(res, bb);
}

void
zsexpmods(
	long a,
	long e,
	verylong n,
	verylong *bb
	)
{
	verylong le = 0;

	zintoz(e, &le);
	zexpmods(a, le, n, bb);
}

void
zexpmods(
	long a,
	verylong e,
	verylong n,
	verylong *bb
	)
{
	long i;
	long j;
	long k = 0;
	verylong res = 0;

	if (a < 0 || zscompare(n, a) <= 0 || zscompare(n, 1) <= 0)
		zhalt("zexpmods: bad args");

	if (!e)
	{
		zone(bb);
		return;
	}

	if ((i = e[0]) < 0)
		zhalt("negative exponent in zexpmods");

	if ((i == 1) && (!(e[1]))) {
		zone(bb);
		return;
	}

	if (!a)
	{
		zzero(bb);
		return;
	}

	zintoz(a, &res);
	for (; i; i--)
	{
		j = e[i];
		if (!k)
		{
			k = 1;
			while ((k << 1) <= j)
				k <<= 1;
		}
		while (k >>= 1)
		{
			zsqmod(res, n, &res);
			if (j & k)
				zsmulmod(res, a, n, &res);
		}
		k = NTL_RADIX;
	}
	zcopy(res, bb);
}



void
z2mul(
	verylong n,
	verylong *rres
	)
{
	long sn;
	long i;
	long carry = 0;
	verylong res = *rres;

	if (!n)
	{
		zzero(rres);
		return;
	}
	if ((!n[1]) && (n[0] == 1))
	{
		zzero(rres);
		return;
	}
	if ((sn = n[0]) < 0)
		sn = -sn;
	zsetlength(&res, sn + 1);
	if (n == *rres) n = res;
	*rres = res;
	for (i = 1; i <= sn; i++)
	{
		if ((res[i] = (n[i] << 1) + carry) >= NTL_RADIX)
		{
			res[i] -= NTL_RADIX;
			carry = 1;
		}
		else
			carry = 0;
	}
	if (carry)
		res[++sn] = 1;
	if (n[0] < 0)
		res[0] = -sn;
	else
		res[0] = sn;
}

long 
z2div(
	verylong n,
	verylong *rres
	)
{
	long sn;
	long i;
	long result;
	verylong res = *rres;

	if ((!n) || ((!n[1]) && (n[0] == 1)))
	{
		zzero(rres);
		return (0);
	}
	if ((sn = n[0]) < 0)
		sn = -sn;
	zsetlength(&res, sn);
	if (n == *rres) n = res;
	*rres = res;
	result = n[1] & 1;
	for (i = 1; i < sn; i++)
	{
		res[i] = (n[i] >> 1);
		if (n[i + 1] & 1)
			res[i] += (NTL_RADIX >> 1);
	}
	if (res[sn] = (n[sn] >> 1))
		res[0] = n[0];
	else if (sn == 1)
	{
		res[0] = 1;
	}
	else if (n[0] < 0)
		res[0] = -sn + 1;
	else
		res[0] = sn - 1;
	return (result);
}

void
zlshift(
	verylong n,
	long k,
	verylong *rres
	)
{
	long big;
	long small;
	long sn;
	long i;
	long cosmall;
	verylong res = *rres;

	if (!n)
	{
		zzero(rres);
		return;
	}
	if ((!n[1]) && (n[0] == 1))
	{
		zzero(rres);
		return;
	}
	if (!k)
	{
		if (n != *rres)
			zcopy(n, rres);
		return;
	}
	if (k < 0)
	{
		if (k < -NTL_MAX_LONG) zhalt("overflow in zlshift");
		zrshift(n, -k, rres);
		return;
	}
	if (k == 1)
	{
		z2mul(n, rres);
		return;
	}
	if ((sn = n[0]) < 0)
		sn = -sn;
	i = sn + (big = k / NTL_NBITS);
	if (small = k - big * NTL_NBITS)
	{
		zsetlength(&res, i + 1);
		if (n == *rres) n = res;
		*rres = res;
		res[i + 1] = n[sn] >> (cosmall = NTL_NBITS - small);
		for (i = sn; i > 1; i--)
			res[i + big] = ((n[i] << small) & NTL_RADIXM) + (n[i - 1] >> cosmall);
		res[big + 1] = (n[1] << small) & NTL_RADIXM;
		for (i = big; i; i--)
			res[i] = 0;
		if (res[sn + big + 1])
			big++;
	}
	else
	{
		zsetlength(&res, i);
		if (n == *rres) n = res;
		*rres = res;
		for (i = sn; i; i--)
			res[i + big] = n[i];
		for (i = big; i; i--)
			res[i] = 0;
	}
	if (n[0] > 0)
		res[0] = n[0] + big;
	else
		res[0] = n[0] - big;
}

void
zrshift(
	verylong n,
	long k,
	verylong *rres
	)
{
	long big;
	long small;
	long sn;
	long i;
	long cosmall;
	verylong res = *rres;

	if (!n)
	{
		zzero(rres);
		return;
	}
	if ((!n[1]) && (n[0] == 1))
	{
		zzero(rres);
		return;
	}
	if (!k)
	{
		if (n != *rres)
			zcopy(n, rres);
		return;
	}
	if (k < 0)
	{
		if (k < -NTL_MAX_LONG) zhalt("overflow in zrshift");
		zlshift(n, -k, rres);
		return;
	}
	if (k == 1)
	{
		z2div(n, rres);
		return;
	}
	big = k / NTL_NBITS;
	small = k - big * NTL_NBITS;
	if ((sn = n[0]) < 0)
		sn = -sn;
	if ((big >= sn) || 
            ((big == sn - 1) && small && (!(n[sn] >> small))))
        /* The microsoft optimizer generates bad code without
           the above test for small != 0 */
	{
		zzero(rres);
		return;
	}
	sn -= big;
	zsetlength(&res, sn);
	if (n == *rres) n = res;
	*rres = res;
	if (small)
	{
		cosmall = NTL_NBITS - small;
		for (i = 1; i < sn; i++)
			res[i] = (n[i + big] >> small) +
				((n[i + big + 1] << cosmall) & NTL_RADIXM);
		if (!(res[sn] = (n[sn + big] >> small)))
			sn--;
	}
	else
		for (i = 1; i <= sn; i++)
			res[i] = n[i + big];
	if (n[0] > 0)
		res[0] = sn;
	else
		res[0] = -sn;
}

long 
zmakeodd(
	verylong *nn
	)
{
	verylong n = *nn;
	long i;
	long shift = 1;

	if (!n || (!n[1] && (n[0] == 1)))
		return (0);
	while (!(n[shift]))
		shift++;
	i = n[shift];
	shift = NTL_NBITS * (shift - 1);
	while (!(i & 1))
	{
		shift++;
		i >>= 1;
	}
	zrshift(n, shift, &n);
	return (shift);
}


long 
zsqrts(
	long n
	)
{
	long a;
	long ndiva;
	long newa;
	verylong ln=0;
	verylong rr=0;

	if (n < 0) 
		zhalt("zsqrts: negative argument");

	if (n <= 0)
		return (0);
	if (n <= 3)
		return (1);
	if (n <= 8)
		return (2);
	if (n >= NTL_RADIX)
	{
		zintoz(n,&ln);
		zsqrt(ln,&rr);
		return(ztoint(rr));
	}
	newa = 3L << (2 * (NTL_NBITSH - 1));
	a = 1L << NTL_NBITSH;
	while (!(n & newa))
	{
		newa >>= 2;
		a >>= 1;
	}
	while (1)
	{
		newa = ((ndiva = n / a) + a) / 2;
		if (newa - ndiva <= 1)
		{
			if (newa * newa <= n)
				return (newa);
			else
				return (ndiva);
		}
		a = newa;
	}
}

void zsqrt(verylong n, verylong *rr)
{
	verylong a = 0;
	verylong ndiva = 0;
	verylong diff = 0;
	verylong r = 0;
	long i;

	if (!n) {
		zzero(rr);
		return;
	}

	if ((i = n[0]) < 0)
		zhalt("negative argument in zsqrt");

	if (i == 1) {
		zintoz(zsqrts(n[1]), rr);
		return;
	}

	zsetlength(&a, i);
	zsetlength(&ndiva, i);
	zsetlength(&diff, i);

	a[(a[0] = (i + 1) / 2)] = zsqrts(n[i]) + 1;
	if (!(i & 1))
		a[a[0]] <<= NTL_NBITSH;

	if (a[a[0]] & NTL_RADIX) {
		a[a[0]] = 0;
		a[0]++;
		a[a[0]] = 1;
	}

	for (i = a[0] - 1; i; i--)
		a[i] = 0;

	while (1) {
		zdiv(n, a, &ndiva, &r);
		zadd(a, ndiva, &r);
		zrshift(r, 1, &r);
		if (zcompare(r, ndiva) <= 0) 
			goto done;

		zsubpos(r, ndiva, &diff);
		if ((diff[0] == 1) && (diff[1] <= 1)) {
			zsq(r, &diff);
			if (zcompare(diff, n) > 0)
				zcopy(ndiva, &r);

			goto done;
		}
		zcopy(r, &a);
	}
done:
	zcopy(r, rr);
}



void
zgcd(
	verylong mm1,
	verylong mm2,
	verylong *rres
	)
{
	long agrb;
	long shibl;
	verylong aa = 0;
	verylong bb = 0;
	verylong cc = 0;
	verylong a;
	verylong b;
	verylong c;
	verylong d;
	long m1negative;
	long m2negative;

	/* ziszero is necessary here and below to fix an
	   an aliasing bug in LIP */

	if (ziszero(mm1))
	{
		if (mm2 != *rres)
			zcopy(mm2,rres);
		zabs(rres);
		return;
	}

	if (ziszero(mm2))
	{
		if (mm1 != *rres)
			zcopy(mm1,rres);
		zabs(rres);
		return;
	}

	if (mm1 == mm2)
	{
		if (mm1 != *rres)
			zcopy(mm1, rres);
		zabs(rres);
		return;
	}

	if (m1negative = (mm1[0] < 0))
		mm1[0] = -mm1[0];
	if (m2negative = (mm2[0] < 0))
		mm2[0] = -mm2[0];

	if ((agrb = mm1[0]) < mm2[0])
		agrb = mm2[0];
	zsetlength(&aa, agrb+1); 
	zsetlength(&bb, agrb+1);
	zsetlength(&cc, agrb+1);
	if (mm1[0] != mm2[0])
	{
		if (mm1[0] > mm2[0])
		{
			zcopy(mm2, &aa);
			zmod(mm1, aa, &bb);
		}
		else
		{
			zcopy(mm1, &aa);
			zmod(mm2, aa, &bb);
		}
		if (!(bb[1]) && (bb[0] == 1))
		{
			a = aa;
			goto done;
		}
	}
	else
	{
		zcopy(mm1, &aa);
		zcopy(mm2, &bb);
	}
	if ((agrb = zmakeodd(&aa)) < (shibl = zmakeodd(&bb)))
		shibl = agrb;
	if (!(agrb = zcompare(aa, bb)))
	{
		a = aa;
		goto endshift;
	}
	else if (agrb < 0)
	{
		a = bb;
		b = aa;
	}
	else
	{
		a = aa;
		b = bb;
	}
	c = cc;
	zsubpos(a, b, &c);
	do
	{
		zmakeodd(&c);
		if (!(agrb = zcompare(b, c)))
		{
			a = b;
			goto endshift;
		}
		else if (agrb > 0)
		{
			a = b;
			b = c;
			c = a;
		}
		else
		{
			d = a;
			a = c;
			c = d;
		}
		zsubpos(a, b, &c);
	} while (c[1] || c[0] != 1);
endshift:
	zlshift(a, shibl, &a);
done:
	if (m1negative)
		mm1[0] = -mm1[0];
	if (m2negative)
		mm2[0] = -mm2[0];
	zcopy(a, rres);
}


long zsign(verylong a)
{
	if (!a)
	{
		return (0);
	}
	if (a[0] < 0)
		return (-1);
	if (a[0] > 1)
		return (1);
	if (a[1])
		return (1);
	return (0);
}

void zabs(verylong *pa)
{
	verylong a = *pa;

	if (!a)
		return;
	if (a[0] < 0)
		a[0] = (-a[0]);
}

long 
z2logs(
	long aa
	)
{
	long i = 0;
	unsigned long a;

	if (aa < 0)
		a = -aa;
	else
		a = aa;

	while (a>=256)
		i += 8, a >>= 8;
	if (a >=16)
		i += 4, a >>= 4;
	if (a >= 4)
		i += 2, a >>= 2;
	if (a >= 2)
		i += 2;
	else if (a >= 1)
		i++;
	return (i);
}

long 
z2log(
	verylong a
	)
{
	long la;

	if (!a)
		return (0);
	la = (a[0] > 0 ? a[0] : -a[0]);
	return ( NTL_NBITS * (la - 1) + z2logs(a[la]) );
}



long
zscompare(
	verylong a,
	long b
	)
{
	if (!b) 
		return zsign(a);
	else {
		verylong c = 0;
		zintoz(b, &c);
		return (zcompare(a, c));
	}
}

void
zswap(
	verylong *a,
	verylong *b
	)
{
	verylong c;

	if ((*a && ((*a)[-1] & 1)) || (*b && ((*b)[-1] & 1))) {
		verylong t = 0; 
		zcopy(*a, &t);
		zcopy(*b, a);
		zcopy(t, b);
		return;
	}

	c = *a;
	*a = *b;
	*b = c;
}

long
ziszero(
	verylong a
	)
{
	if (!a) return (1);
	if (a[1]) return (0);
	if (a[0]==1) return (1);
	return (0);
}

long
zodd(
	verylong a
	)
{
	if (!a) return (0);
	return (a[1]&1);
}


long
zbit(
	verylong a,
	long p
	)
{
        long bl;
        long wh;
        long sa;
	if (p < 0 || !a) return 0;
	bl = (p/NTL_NBITS);
        wh = 1L << (p - NTL_NBITS*bl);
	bl ++;
        sa = a[0];
        if (sa < 0) sa = -sa;
        if (sa < bl) return (0);
	if (a[bl] & wh) return (1);
	return (0);
}


void
zlowbits(
	verylong a,
	long b,
        verylong *cc
        )
{
        verylong c = *cc;
	long bl;
	long wh;
	long sa;

	if (ziszero(a) || (b<=0)) {
		zzero(cc);
		return;
	}

        bl = b/NTL_NBITS;
        wh = b - NTL_NBITS*bl;
	if (wh != 0) 
		bl++;
	else
		wh = NTL_NBITS;

	sa = a[0];
	if (sa < 0) sa = -sa;

	if (sa < bl) {
		zcopy(a,cc);
		zabs(cc);
		return;
	}

        zsetlength(&c, bl);
        if (a == *cc) a = c;
        *cc = c;
	for (sa=1; sa<bl; sa++)
		c[sa] = a[sa];
	c[bl] = a[bl]&((1L<<wh)-1);
	while ((bl>1) && (!c[bl]))
		bl --;
	c[0] = bl;
}


long zslowbits(verylong a, long p)
{
   verylong x = 0;

   zlowbits(a, p, &x);

   return ztoint(x);
}




long
zweights(
	long aa
	)
{
	unsigned long a;
	long res = 0;
	if (aa < 0) 
		a = -aa;
	else
		a = aa;
   
	while (a) {
		if (a & 1) res ++;
		a >>= 1;
	}
	return (res);
}

long
zweight(
        verylong a
        )
{
	long i;
	long res = 0;
	if (!a) return (0);
	i = a[0];
	if (i<0) i = -i;
	for (;i;i--)
		res += zweights(a[i]);
	return (res);
}



void
zand(
	verylong a,
	verylong b,
	verylong *cc
	)

{
	verylong c = *cc;
	long sa;
	long sb;
	long sm;
	if (ziszero(a) || ziszero(b)) {
		zzero(cc);
		return;
	}
	sa = a[0];
	if (sa < 0) sa = -sa;
	sb = b[0];
	if (sb < 0) sb = -sb;
	sm = (sa > sb ? sb : sa );
	zsetlength(&c, sm);
	if (a == *cc) a = c;
	if (b == *cc) b = c;
	*cc = c;
	for (sa = 1; sa <= sm; sa ++) 
		c[sa] = a[sa] & b[sa];
	while ((sm > 1) && (!(c[sm])))
		sm --;
	c[0] = sm;
}

void
zxor(
	verylong a,
	verylong b,
	verylong *cc
	)
{
        verylong c = *cc;
        long sa;
        long sb;
        long sm;
        long la;
        long i;
	if (ziszero(a)) {
		zcopy(b,cc);
		zabs(cc);
		return;
	}
	if (ziszero(b)) {
		zcopy(a,cc);
		zabs(cc);
		return;
	}
        sa = a[0];
        if (sa < 0) sa = -sa;
        sb = b[0];
        if (sb < 0) sb = -sb;
	if (sa > sb) {
		la = sa;
		sm = sb;
	} else {
		la = sb;
		sm = sa;
	}
        zsetlength(&c, la);
        if (a == *cc) a = c;
        if (b == *cc) b = c;
        *cc = c;
        for (i = 1; i <= sm; i ++)
                c[i] = a[i] ^ b[i];
	if (sa > sb)
		for (;i <= la; i++) c[i] = a[i];
	else
		for (;i <= la; i++) c[i] = b[i];
        while ((la > 1) && (!(c[la])))
                la --;
        c[0] = la;
}

void
zor(
	verylong a,
	verylong b,
        verylong *cc
        )
{
        verylong c = *cc;
        long sa;
        long sb;
        long sm;
        long la;
        long i;
	if (ziszero(a)) {
		zcopy(b,cc);
		zabs(cc);
		return;
	}
	if (ziszero(b)) {
		zcopy(a,cc);
		zabs(cc);
		return;
	}
        sa = a[0];
        if (sa < 0) sa = -sa;
        sb = b[0];
        if (sb < 0) sb = -sb;
        if (sa > sb) {
                la = sa;
                sm = sb;
        } else {
                la = sb;
                sm = sa;
        }
        zsetlength(&c, la);
        if (a == *cc) a = c;
        if (b == *cc) b = c;
        *cc = c;
        for (i = 1; i <= sm; i ++)
                c[i] = a[i] | b[i];
        if (sa > sb)
                for (;i <= la; i++) c[i] = a[i];
        else
                for (;i <= la; i++) c[i] = b[i];
        c[0] = la;
}

long
zsetbit(
        verylong *a,
        long b
        )
{
	long bl;
	long wh;
	long sa;

	if (b<0) zhalt("zsetbit: negative index");

	if (ziszero(*a)) {
		zintoz(1,a);
		zlshift(*a,b,a);
		return (0);
	}

	bl = (b/NTL_NBITS);
	wh = 1L << (b - NTL_NBITS*bl);
	bl ++;
	sa = (*a)[0];
	if (sa<0) sa = -sa;
	if (sa >= bl) {
		sa = (*a)[bl] & wh;
		(*a)[bl] |= wh;
		if (sa) return (1);
		return (0);
	} else {
		zsetlength(a,bl);
		sa ++;
		for (;sa<=bl;sa++) (*a)[sa]=0;
		if ((*a)[0] < 0)
			(*a)[0] = -bl;
		else (*a)[0] = bl;
		(*a)[bl] |= wh;
		return (0);
	}
}

long
zswitchbit(
        verylong *a,
        long p
        )
{
        long bl;
        long wh;
        long sa;

        if (p < 0) zhalt("zswitchbit: negative index");

        if (ziszero(*a)) {
		zintoz(1,a);
		zlshift(*a,p,a);
		return (0);
	}

        bl = (p/NTL_NBITS);
        wh = 1L << (p - NTL_NBITS*bl);
        bl ++;
        sa = (*a)[0];
        if (sa < 0) sa = -sa;
        if ((sa < bl) || (!((*a)[bl] & wh))) {
		zsetbit(a,p);
		return (0);
	}
	(*a)[bl] ^= wh;
        while ((sa>1) && (!(*a)[sa]))
		sa --;
	if ((*a)[0] > 0) (*a)[0] = sa;
	else (*a)[0] = -sa;
        return (1);
}

