

#ifndef NTL_ZZ__H
#define NTL_ZZ__H



/********************************************************

   LIP INTERFACE 

   The class ZZ implements signed, arbitrary length integers.

**********************************************************/


#include <NTL/lip.h>
#include <NTL/tools.h>


class ZZ {
public:

long *rep; // This is currently public for "emergency" situations
           // May be private in future versions.


ZZ() 
// initial value is 0.

{ rep = 0; }


ZZ(INIT_SIZE_TYPE, long k)
// initial value is 0, but space is pre-allocated so that numbers
// x with x.size() <= k can be stored without re-allocation.
// Call with ZZ(INIT_SIZE, k).
// The purpose for the INIT_SIZE argument is to prevent automatic
// type conversion from long to ZZ, which would be tempting, but wrong.


{
   rep = 0;
   zsetlength(&rep, k); 
}

ZZ(const ZZ& a)
// initial value is a.

{
   rep = 0;
   zcopy(a.rep, &rep);
}


ZZ(INIT_VAL_TYPE, long a) { rep = 0; zintoz(a, &rep); }
ZZ(INIT_VAL_TYPE, int a) { rep = 0; zintoz(a, &rep); }

inline ZZ(INIT_VAL_TYPE, const char *);
inline ZZ(INIT_VAL_TYPE, float);
inline ZZ(INIT_VAL_TYPE, double);


ZZ& operator=(const ZZ& a) { zcopy(a.rep, &rep); return *this; }

ZZ& operator=(long a) { zintoz(a, &rep); return *this; }


~ZZ() { zfree(&rep); }

void kill()
// force the space held by this ZZ to be released.
// The value then becomes 0.

{ zfree(&rep); }

void SetSize(long k)
// pre-allocates space for k-digit numbers (base ZZ_RADIX);  
// does not change the value.

{ zsetlength(&rep, k); }

long size() const; 
// returns the number of (ZZ_NBIT-bit) digits of |a|; the size of 0 is 0.

static const ZZ& zero();


ZZ(ZZ& x, INIT_TRANS_TYPE) { rep = x.rep; x.rep = 0; }
// used to cheaply hand off memory management of return value,
// without copying, assuming compiler implements the
// "return value optimization"

};



const ZZ& ZZ_expo(long e);


inline void clear(ZZ& x)
// x = 0

   { zzero(&x.rep); }

inline void set(ZZ& x)
// x = 1

   { zsetlength(&x.rep, 1); x.rep[0] = 1; x.rep[1] = 1; }


inline void swap(ZZ& x, ZZ& y)
// swap the values of x and y (swaps pointers only)

   { zswap(&x.rep, &y.rep); }


double log(const ZZ& a);




/**********************************************************

   Conversion routines.

***********************************************************/



inline void conv(ZZ& x, const ZZ& a) { x = a; }
inline ZZ to_ZZ(const ZZ& a) { return a; }

inline void conv(ZZ& x, long a) { zintoz(a, &x.rep); }
inline ZZ to_ZZ(long a) { return ZZ(INIT_VAL, a); }


inline void conv(ZZ& x, int a) { zintoz(long(a), &x.rep); }
inline ZZ to_ZZ(int a) { return ZZ(INIT_VAL, a); }

void conv(ZZ& x, const char *s);
inline ZZ::ZZ(INIT_VAL_TYPE, const char *s) {  rep = 0; conv(*this, s); }
inline ZZ to_ZZ(const char *s) { return ZZ(INIT_VAL, s); }

inline void conv(ZZ& x, double a) { zdoubtoz(a, &x.rep); }
inline ZZ::ZZ(INIT_VAL_TYPE, double a) { rep = 0; conv(*this, a); }
inline ZZ to_ZZ(double a) { return ZZ(INIT_VAL, a); }

inline void conv(ZZ& x, float a) { zdoubtoz(double(a), &x.rep); }
inline ZZ::ZZ(INIT_VAL_TYPE, float a) { rep = 0; conv(*this, a); }
inline ZZ to_ZZ(float a) { return ZZ(INIT_VAL, a); }

inline void conv(long& x, const ZZ& a) { x = ztoint(a.rep); }
inline long to_long(const ZZ& a)  { return ztoint(a.rep); }

inline void conv(int& x, const ZZ& a) { x = int(ztoint(a.rep)); }
inline int to_int(const ZZ& a)  { return int(ztoint(a.rep)); }

inline void conv(double& x, const ZZ& a) { x = zdoub(a.rep); }
inline double to_double(const ZZ& a) { return zdoub(a.rep); }

inline void conv(float& x, const ZZ& a) { x = float(zdoub(a.rep)); }
inline float to_float(const ZZ& a) { return float(zdoub(a.rep)); }

void ZZFromBytes(ZZ& x, const unsigned char *p, long n);
inline ZZ ZZFromBytes(const unsigned char *p, long n)
   { ZZ x; ZZFromBytes(x, p, n); NTL_OPT_RETURN(ZZ, x); }

void BytesFromZZ(unsigned char *p, const ZZ& a, long n);




// ****** comparisons


inline long sign(const ZZ& a)
// returns the sign of a (-1, 0, or 1).

   { return zsign(a.rep); }


inline long compare(const ZZ& a, const ZZ& b)
// returns the sign of a-b (-1, 0, or 1).

{
   return zcompare(a.rep, b.rep);
}

inline long IsZero(const ZZ& a)
// zero test

   { return ziszero(a.rep); }


long IsOne(const ZZ& a);
// test for 1
   

/* the usual comparison operators */

inline long operator==(const ZZ& a, const ZZ& b)
  { return zcompare(a.rep, b.rep) == 0; }
inline long operator!=(const ZZ& a, const ZZ& b)
  { return zcompare(a.rep, b.rep) != 0; }
inline long operator<(const ZZ& a, const ZZ& b)
  { return zcompare(a.rep, b.rep) < 0; }
inline long operator>(const ZZ& a, const ZZ& b)
  { return zcompare(a.rep, b.rep) > 0; }
inline long operator<=(const ZZ& a, const ZZ& b)
  { return zcompare(a.rep, b.rep) <= 0; }
inline long operator>=(const ZZ& a, const ZZ& b)
  { return zcompare(a.rep, b.rep) >= 0; }

/* single-precision versions of the above */

inline long compare(const ZZ& a, long b) { return zscompare(a.rep, b); }
inline long compare(long a, const ZZ& b) { return -zscompare(b.rep, a); }

inline long operator==(const ZZ& a, long b) { return zscompare(a.rep, b) == 0; }
inline long operator!=(const ZZ& a, long b) { return zscompare(a.rep, b) != 0; }
inline long operator<(const ZZ& a, long b) { return zscompare(a.rep, b) < 0; }
inline long operator>(const ZZ& a, long b) { return zscompare(a.rep, b) > 0; }
inline long operator<=(const ZZ& a, long b) { return zscompare(a.rep, b) <= 0; }
inline long operator>=(const ZZ& a, long b) { return zscompare(a.rep, b) >= 0; }


inline long operator==(long a, const ZZ& b) { return b == a; }
inline long operator!=(long a, const ZZ& b) { return b != a; }
inline long operator<(long a, const ZZ& b) { return b > a; }
inline long operator>(long a, const ZZ& b) { return b < a; }
inline long operator<=(long a, const ZZ& b) { return b >= a; }
inline long operator>=(long a, const ZZ& b) { return b <= a; }

/**************************************************

                 Addition

**************************************************/


inline void add(ZZ& x, const ZZ& a, const ZZ& b)
// x = a + b

   { zadd(a.rep, b.rep, &x.rep); }

inline void sub(ZZ& x, const ZZ& a, const ZZ& b)
// x = a - b

   { zsub(a.rep, b.rep, &x.rep); }

inline void SubPos(ZZ& x, const ZZ& a, const ZZ& b)
// x = a - b;  assumes a >= b >= 0.

   { zsubpos(a.rep, b.rep, &x.rep); }

inline void negate(ZZ& x, const ZZ& a)
// x = -a

   { zcopy(a.rep, &x.rep); znegate(&x.rep); }

inline void abs(ZZ& x, const ZZ& a)
// x = |a|
{ zcopy(a.rep, &x.rep); zabs(&x.rep); }


/* single-precision versions of the above */

inline void add(ZZ& x, const ZZ& a, long b)
   { zsadd(a.rep, b, &x.rep); }

inline void add(ZZ& x, long a, const ZZ& b) { add(x, b, a); }


void sub(ZZ& x, const ZZ& a, long b);
void sub(ZZ& x, long a, const ZZ& b);

/* operator/function notation */

inline ZZ operator+(const ZZ& a, const ZZ& b) 
  { ZZ x; add(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator+(const ZZ& a, long b) 
  { ZZ x; add(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator+(long  a, const ZZ& b) 
  { ZZ x; add(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator-(const ZZ& a, const ZZ& b) 
  { ZZ x; sub(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator-(const ZZ& a, long b) 
  { ZZ x; sub(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator-(long  a, const ZZ& b) 
  { ZZ x; sub(x, a, b); NTL_OPT_RETURN(ZZ, x); } 

inline ZZ operator-(const ZZ& a)
  { ZZ x; negate(x, a); NTL_OPT_RETURN(ZZ, x); }

inline ZZ abs(const ZZ& a)
  { ZZ x; abs(x, a); NTL_OPT_RETURN(ZZ, x); }

/* op= notation */

inline ZZ& operator+=(ZZ& x, const ZZ& a)
  { add(x, x, a); return x; }

inline ZZ& operator+=(ZZ& x, long a)
  { add(x, x, a); return x; }

inline ZZ& operator-=(ZZ& x, const ZZ& a)
  { sub(x, x, a); return x; }

inline ZZ& operator-=(ZZ& x, long a)
  { sub(x, x, a); return x; }

/* inc/dec */

inline ZZ& operator++(ZZ& x) { add(x, x, 1); return x; }

inline void operator++(ZZ& x, int) { add(x, x, 1); }

inline ZZ& operator--(ZZ& x) { add(x, x, -1); return x; }

inline void operator--(ZZ& x, int) { add(x, x, -1); }



/*******************************************************

                 Multiplication.

********************************************************/

inline void mul(ZZ& x, const ZZ& a, const ZZ& b)
// x = a * b

   { zmul(a.rep, b.rep, &x.rep); }


inline void sqr(ZZ& x, const ZZ& a)
// x = a*a

   { zsq(a.rep, &x.rep); }

inline ZZ sqr(const ZZ& a)
   { ZZ x; sqr(x, a); NTL_OPT_RETURN(ZZ, x); }


/* single-precision versions */

inline void mul(ZZ& x, const ZZ& a, long b)
   { zsmul(a.rep, b, &x.rep); }

inline void mul(ZZ& x, long a, const ZZ& b)
    { mul(x, b, a); }

/* operator notation */

inline ZZ operator*(const ZZ& a, const ZZ& b)
  { ZZ x; mul(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator*(const ZZ& a, long b)
  { ZZ x; mul(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator*(long a, const ZZ& b)
  { ZZ x; mul(x, a, b); NTL_OPT_RETURN(ZZ, x); }

/* op= notation */

inline ZZ& operator*=(ZZ& x, const ZZ& a)
  { mul(x, x, a); return x; }

inline ZZ& operator*=(ZZ& x, long a)
  { mul(x, x, a); return x; }


void MultiMul(ZZ& x, long n, const ZZ* a, const long* b, long size);
// x = sum_{i=0}^{n-1} a[i]*b[i].
// It is assumed that a is nonnegative, and that
// 0 <= b[i] < 2^{NTL_NBITS} for each i.
// ALIAS RESTRICTION: inputs can not alias outputs.
// The parameter size should be computed as an upper bound
// on the size of the result.
// This routine is heavily used in chinese remaindering,
// and therefore must be as fast as possible.




/*******************************************************

                    Division

*******************************************************/


inline void DivRem(ZZ& q, ZZ& r, const ZZ& a, const ZZ& b)
// q = [a/b], r = a - b*q
// |r| < |b|, and if r != 0, sign(r) = sign(b)

   { zdiv(a.rep, b.rep, &q.rep, &r.rep); }



inline void div(ZZ& q, const ZZ& a, const ZZ& b)
// q = a/b

   { zdiv(a.rep, b.rep, &q.rep, 0); }

inline void rem(ZZ& r, const ZZ& a, const ZZ& b)
// r = a%b

   { zmod(a.rep, b.rep, &r.rep); }


inline void QuickRem(ZZ& r, const ZZ& b)
// r = r%b
// assumes b > 0 and r >=0
// division is performed in place and may cause r to be re-allocated.

   { zquickmod(&r.rep, b.rep); }

long divide(ZZ& q, const ZZ& a, const ZZ& b);
// if b | a, sets q = a/b and returns 1; otherwise returns 0.

long divide(const ZZ& a, const ZZ& b);
// if b | a, returns 1; otherwise returns 0.


/* non-standard single-precision versions */

inline long DivRem(ZZ& q, const ZZ& a, long b)
   { return zsdiv(a.rep, b, &q.rep); } 

inline long rem(const ZZ& a, long b)
   { return zsmod(a.rep, b); }


/* single precision versions */

inline void div(ZZ& q, const ZZ& a, long b)
   { (void) zsdiv(a.rep, b, &q.rep); }


long divide(ZZ& q, const ZZ& a, long b);
// if b | a, sets q = a/b and returns 1; otherwise returns 0.

long divide(const ZZ& a, long b);
// if b | a, returns 1; otherwise returns 0.


#if (defined(NTL_SINGLE_MUL))

inline void MultiRem(long* r, long n, const ZZ& a, 
                     const long* d, double **tbl)

   { zmultirem2(a.rep, n, (long *) d, tbl, r); }

#elif (defined(NTL_TBL_REM))

inline void MultiRem(long* r, long n, const ZZ& a, 
                     const long* d, long **tbl)

   { zmultirem3(a.rep, n, (long *) d, tbl, r); }

#else 

inline void MultiRem(long* r, long n, const ZZ& a, const long* d)

   { zmultirem(a.rep, n, (long *) d, r); }

#endif

// computes r[i] = a % d[i], i = 0..n-1.
// assumes a >= 0 and 0 < d[i] < 2^{NTL_NBITS}
// This must be fast, as it is critical in chinese remaindering.


inline ZZ operator/(const ZZ& a, const ZZ& b)
   { ZZ x; div(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator/(const ZZ& a, long b)
   { ZZ x; div(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator%(const ZZ& a, const ZZ& b)
   { ZZ x; rem(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline long operator%(const ZZ& a, long b)
   { return rem(a, b); }

inline ZZ& operator/=(ZZ& x, const ZZ& b)
   { div(x, x, b); return x; } 

inline ZZ& operator/=(ZZ& x, long b)
   { div(x, x, b); return x; } 

inline ZZ& operator%=(ZZ& x, const ZZ& b)
   { rem(x, x, b); return x; } 


/**********************************************************

                        GCD's

***********************************************************/


inline void GCD(ZZ& d, const ZZ& a, const ZZ& b)
// d = gcd(a, b)

   { zgcd(a.rep, b.rep, &d.rep); }

inline ZZ GCD(const ZZ& a, const ZZ& b)
   { ZZ x; GCD(x, a, b); NTL_OPT_RETURN(ZZ, x); }


inline void XGCD(ZZ& d, ZZ& s, ZZ& t, const ZZ& a, const ZZ& b)
//  d = gcd(a, b) = a*s + b*t;

   { zexteucl(a.rep, &s.rep, b.rep, &t.rep, &d.rep); }

// single-precision versions
long GCD(long a, long b);

void XGCD(long& d, long& s, long& t, long a, long b);







/************************************************************

                      Bit Operations

*************************************************************/


inline void LeftShift(ZZ& x, const ZZ& a, long k)
// x = (a << k), k < 0 => RightShift

   { zlshift(a.rep, k, &x.rep); }

inline ZZ LeftShift(const ZZ& a, long k)
   { ZZ x; LeftShift(x, a, k); NTL_OPT_RETURN(ZZ, x); }


inline void RightShift(ZZ& x, const ZZ& a, long k)
// x = (a >> k), k < 0 => LeftShift

   { zrshift(a.rep, k, &x.rep); }

inline ZZ RightShift(const ZZ& a, long k)
   { ZZ x; RightShift(x, a, k); NTL_OPT_RETURN(ZZ, x); }

#ifndef NTL_TRANSITION

inline ZZ operator>>(const ZZ& a, long n)
   { ZZ x; RightShift(x, a, n); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator<<(const ZZ& a, long n)
   { ZZ x; LeftShift(x, a, n); NTL_OPT_RETURN(ZZ, x); }

inline ZZ& operator<<=(ZZ& x, long n)
   { LeftShift(x, x, n); return x; }

inline ZZ& operator>>=(ZZ& x, long n)
   { RightShift(x, x, n); return x; }

#endif


inline long MakeOdd(ZZ& x)
// removes factors of 2 from x, returns the number of 2's removed
// returns 0 if x == 0

   { return zmakeodd(&x.rep); }


inline long IsOdd(const ZZ& a)
// returns 1 if a is odd, otherwise 0

   { return zodd(a.rep); }


inline long NumBits(const ZZ& a)
// returns the number of bits in |a|; NumBits(0) = 0
   { return z2log(a.rep); }



inline long bit(const ZZ& a, long k)
// returns bit k of a, 0 being the low-order bit

   { return zbit(a.rep, k); }

long digit(const ZZ& a, long k);
// returns k-th digit of |a|, 0 being the low-order digit.

inline void trunc(ZZ& x, const ZZ& a, long k)
// puts k low order bits of |a| into x

   { zlowbits(a.rep, k, &x.rep); }

inline ZZ trunc_ZZ(const ZZ& a, long k)
   { ZZ x; trunc(x, a, k); NTL_OPT_RETURN(ZZ, x); }

inline long trunc_long(const ZZ& a, long k)
// returns k low order bits of |a|

   { return zslowbits(a.rep, k); }

inline long SetBit(ZZ& x, long p)
// returns original value of p-th bit of |a|, and replaces
// p-th bit of a by 1 if it was zero;
// error if p < 0 */

   { return zsetbit(&x.rep, p); }

inline long SwitchBit(ZZ& x, long p)
// returns original value of p-th bit of |a|, and switches
// the value of p-th bit of a;
// p starts counting at 0;
//   error if p < 0 */

   { return zswitchbit(&x.rep, p); }

inline long weight(long a)
// returns Hamming weight of |a|

   { return zweights(a); }

inline long weight(const ZZ& a)
// returns Hamming weight of |a|

   { return zweight(a.rep); }

inline void bit_and(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| AND |b|

   { zand(a.rep, b.rep, &x.rep); }

void bit_and(ZZ& x, const ZZ& a, long b);
inline void bit_and(ZZ& x, long a, const ZZ& b)
   { bit_and(x, b, a); }


inline void bit_or(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| OR |b|

   { zor(a.rep, b.rep, &x.rep); }

void bit_or(ZZ& x, const ZZ& a, long b);
inline void bit_or(ZZ& x, long a, const ZZ& b)
   { bit_or(x, b, a); }

inline void bit_xor(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| XOR |b|

   { zxor(a.rep, b.rep, &x.rep); }

void bit_xor(ZZ& x, const ZZ& a, long b);
inline void bit_xor(ZZ& x, long a, const ZZ& b)
   { bit_xor(x, b, a); }


inline ZZ operator&(const ZZ& a, const ZZ& b)
   { ZZ x; bit_and(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator&(const ZZ& a, long b)
   { ZZ x; bit_and(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator&(long a, const ZZ& b)
   { ZZ x; bit_and(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator|(const ZZ& a, const ZZ& b)
   { ZZ x; bit_or(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator|(const ZZ& a, long b)
   { ZZ x; bit_or(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator|(long a, const ZZ& b)
   { ZZ x; bit_or(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator^(const ZZ& a, const ZZ& b)
   { ZZ x; bit_xor(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator^(const ZZ& a, long b)
   { ZZ x; bit_xor(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ operator^(long a, const ZZ& b)
   { ZZ x; bit_xor(x, a, b); NTL_OPT_RETURN(ZZ, x); }

inline ZZ& operator&=(ZZ& x, const ZZ& b) 
   { bit_and(x, x, b); return x; }

inline ZZ& operator&=(ZZ& x, long b) 
   { bit_and(x, x, b); return x; }

inline ZZ& operator|=(ZZ& x, const ZZ& b) 
   { bit_or(x, x, b); return x; }

inline ZZ& operator|=(ZZ& x, long b) 
   { bit_or(x, x, b); return x; }

inline ZZ& operator^=(ZZ& x, const ZZ& b) 
   { bit_xor(x, x, b); return x; }

inline ZZ& operator^=(ZZ& x, long b) 
   { bit_xor(x, x, b); return x; }



long NumBits(long a);

long bit(long a, long k);

long NextPowerOfTwo(long m);
// returns least nonnegative k such that 2^k >= m

inline
long NumBytes(const ZZ& a)
   { return (NumBits(a)+7)/8; }

inline
long NumBytes(long a)
   { return (NumBits(a)+7)/8; }


/***********************************************************

                  Psuedo-random Numbers

************************************************************/


void SetSeed(const ZZ& s);
// initialize random number generator

void RandomBnd(ZZ& x, const ZZ& n);
// x = "random number" in the range 0..n-1, or 0  if n <= 0

inline ZZ RandomBnd(const ZZ& n)
   { ZZ x; RandomBnd(x, n); NTL_OPT_RETURN(ZZ, x); }


void RandomLen(ZZ& x, long NumBits);
// x = "random number" with precisely NumBits bits.


inline ZZ RandomLen_ZZ(long NumBits)
   { ZZ x; RandomLen(x, NumBits); NTL_OPT_RETURN(ZZ, x); }


void RandomBits(ZZ& x, long NumBits);
// x = "random number", 0 <= x < 2^NumBits 

inline ZZ RandomBits_ZZ(long NumBits)
   { ZZ x; RandomBits(x, NumBits); NTL_OPT_RETURN(ZZ, x); }


// single-precision version of the above

long RandomBnd(long n);

long RandomLen_long(long l);

long RandomBits_long(long l);

unsigned long RandomWord();



/**********************************************************

             Incremental Chinese Remaindering

***********************************************************/

long CRT(ZZ& a, ZZ& p, const ZZ& A, const ZZ& P);
long CRT(ZZ& a, ZZ& p, long A, long P);
// 0 <= A < P, (p, P) = 1;
// computes b such that b = a mod p, b = A mod p,
//   and -p*P/2 < b <= p*P/2;
// sets a = b, p = p*P, and returns 1 if a's value
//   has changed, otherwise 0

long CRTInRange(const ZZ& gg, const ZZ& aa);
// an auxilliary routine used by newer CRT routines to maintain
// backward compatability.


/**********************************************************

                  Rational Reconstruction

***********************************************************/

inline
long ReconstructRational(ZZ& a, ZZ& b, const ZZ& u, const ZZ& m, 
                         const ZZ& a_bound, ZZ& b_bound)
{
   return zxxratrecon(u.rep, m.rep, a_bound.rep, b_bound.rep, &a.rep, &b.rep);

}




/************************************************************

                      Primality Testing 

*************************************************************/


void GenPrime(ZZ& n, long l, long err = 80);
inline ZZ GenPrime_ZZ(long l, long err = 80) 
{ ZZ x; GenPrime(x, l, err); NTL_OPT_RETURN(ZZ, x); }
long GenPrime_long(long l, long err = 80);
// This generates a random prime n of length l so that the
// probability of erroneously returning a composite is bounded by 2^(-err).

long ProbPrime(const ZZ& n, long NumTrials = 10);
// tests if n is prime;  performs a little trial division,
// followed by a single-precision MillerWitness test, followed by
// up to NumTrials general MillerWitness tests.

long MillerWitness(const ZZ& n, const ZZ& w);
// Tests if w is a witness to primality a la Miller.
// Assumption: n is odd and positive, 0 <= w < n.

void RandomPrime(ZZ& n, long l, long NumTrials=10);
// n =  random l-bit prime

inline ZZ RandomPrime_ZZ(long l, long NumTrials=10)
   { ZZ x; RandomPrime(x, l, NumTrials); NTL_OPT_RETURN(ZZ, x); }

void NextPrime(ZZ& n, const ZZ& m, long NumTrials=10);
// n = smallest prime >= m.

inline ZZ NextPrime(const ZZ& m, long NumTrials=10)
   { ZZ x; NextPrime(x, m, NumTrials); NTL_OPT_RETURN(ZZ, x); }

// single-precision versions

long ProbPrime(long n, long NumTrials = 10);


long RandomPrime_long(long l, long NumTrials=10);

long NextPrime(long l, long NumTrials=10);


/************************************************************

                      Exponentiation

*************************************************************/

void power(ZZ& x, const ZZ& a, long e);
// x = a^e, e >= 0

inline ZZ power(const ZZ& a, long e)
   { ZZ x; power(x, a, e); NTL_OPT_RETURN(ZZ, x); }

inline void power(ZZ& x, long a, long e)
   {  zexps(a, e, &x.rep); }

inline ZZ power_ZZ(long a, long e)
   { ZZ x; power(x, a, e); NTL_OPT_RETURN(ZZ, x); }

long power_long(long a, long e); 

void power2(ZZ& x, long e);

inline ZZ power2_ZZ(long e)
   { ZZ x; power2(x, e); NTL_OPT_RETURN(ZZ, x); }





/*************************************************************

                       Square Roots

**************************************************************/




inline void SqrRoot(ZZ& x, const ZZ& a)
// x = [a^{1/2}], a >= 0

{
   zsqrt(a.rep, &x.rep);
}

inline ZZ SqrRoot(const ZZ& a)
   { ZZ x; SqrRoot(x, a); NTL_OPT_RETURN(ZZ, x); }


inline long SqrRoot(long a) { return zsqrts(a); }
// single-precision version



/***************************************************************

                      Modular Arithmetic

***************************************************************/

// The following routines perform arithmetic mod n, n positive.
// All args (other than exponents) are assumed to be in the range 0..n-1.


void AddMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n);
// x = (a+b)%n

inline ZZ AddMod(const ZZ& a, const ZZ& b, const ZZ& n)
   { ZZ x; AddMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void SubMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
// x = (a-b)%n

   { zsubmod(a.rep, b.rep, n.rep, &x.rep); }

inline ZZ SubMod(const ZZ& a, const ZZ& b, const ZZ& n)
   { ZZ x; SubMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void NegateMod(ZZ& x, const ZZ& a, const ZZ& n)
// x = -a % n

   { zsubmod(0, a.rep, n.rep, &x.rep); }

inline ZZ NegateMod(const ZZ& a, const ZZ& n)
   { ZZ x; NegateMod(x, a, n); NTL_OPT_RETURN(ZZ, x); }

void AddMod(ZZ& x, const ZZ& a, long b, const ZZ& n);
inline ZZ AddMod(const ZZ& a, long b, const ZZ& n)
   { ZZ x; AddMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void AddMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
   { AddMod(x, b, a, n); }
inline ZZ AddMod(long a, const ZZ& b, const ZZ& n)
   { ZZ x; AddMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

void SubMod(ZZ& x, const ZZ& a, long b, const ZZ& n);
inline ZZ SubMod(const ZZ& a, long b, const ZZ& n)
   { ZZ x; SubMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

void SubMod(ZZ& x, long a, const ZZ& b, const ZZ& n);
inline ZZ SubMod(long a, const ZZ& b, const ZZ& n)
   { ZZ x; SubMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void MulMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
// x = (a*b)%n

   { zmulmod(a.rep, b.rep, n.rep, &x.rep); }

inline ZZ MulMod(const ZZ& a, const ZZ& b, const ZZ& n)
   { ZZ x; MulMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void MulMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
// x = (a*b)%n

   { zsmulmod(a.rep, b, n.rep, &x.rep); }

inline ZZ MulMod(const ZZ& a, long b, const ZZ& n)
   { ZZ x; MulMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }

inline void MulMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
   { MulMod(x, b, a, n); }

inline ZZ MulMod(long a, const ZZ& b, const ZZ& n)
   { ZZ x; MulMod(x, a, b, n); NTL_OPT_RETURN(ZZ, x); }


inline void SqrMod(ZZ& x, const ZZ& a, const ZZ& n)
// x = a^2 % n

   { zsqmod(a.rep, n.rep, &x.rep); }

inline ZZ SqrMod(const ZZ& a, const ZZ& n)
   {  ZZ x; SqrMod(x, a, n); NTL_OPT_RETURN(ZZ, x); }

inline void InvMod(ZZ& x, const ZZ& a, const ZZ& n)
// x = a^{-1} mod n, 0 <= x < n
// error is raised occurs if inverse not defined

   { zinvmod(a.rep, n.rep, &x.rep); }

inline ZZ InvMod(const ZZ& a, const ZZ& n)
   {  ZZ x; InvMod(x, a, n); NTL_OPT_RETURN(ZZ, x); }


inline long InvModStatus(ZZ& x, const ZZ& a, const ZZ& n)
// if gcd(a,b) = 1, then ReturnValue = 0, x = a^{-1} mod n
// otherwise, ReturnValue = 1, x = gcd(a, n)

  { return zinv(a.rep, n.rep, &x.rep); }


void PowerMod(ZZ& x, const ZZ& a, const ZZ& e, const ZZ& n);
inline ZZ PowerMod(const ZZ& a, const ZZ& e, const ZZ& n)
   { ZZ x; PowerMod(x, a, e, n); NTL_OPT_RETURN(ZZ, x); }

inline void PowerMod(ZZ& x, const ZZ& a, long e, const ZZ& n)
   { PowerMod(x, a, ZZ_expo(e), n); }

inline ZZ PowerMod(const ZZ& a, long e, const ZZ& n)
   { ZZ x; PowerMod(x, a, e, n); NTL_OPT_RETURN(ZZ, x); }






/*************************************************************

   Jacobi symbol and modular squre roots

**************************************************************/


long Jacobi(const ZZ& a, const ZZ& n);
//  compute Jacobi symbol of a and n;
//  assumes 0 <= a < n, n odd

void SqrRootMod(ZZ& x, const ZZ& a, const ZZ& n);
//  computes square root of a mod n;
//  assumes n is an odd prime, and that a is a square mod n

inline ZZ SqrRootMod(const ZZ& a, const ZZ& n)
   { ZZ x; SqrRootMod(x, a, n); NTL_OPT_RETURN(ZZ, x); }




/*************************************************************


                    Small Prime Generation


*************************************************************/


// primes are generated in sequence, starting at 2, 
// and up until (2*NTL_PRIME_BND+1)^2, which is less than NTL_RADIX.

#if (NTL_NBITS > 30)
#define NTL_PRIME_BND ((1L << 14) - 1)
#else
#define NTL_PRIME_BND ((1L << (NTL_NBITS/2-1)) - 1)
#endif


class PrimeSeq {


char *movesieve;
char *movesieve_mem;
long pindex;
long pshift;
long exhausted;

public:

PrimeSeq();
~PrimeSeq();

long next();
// returns next prime in the sequence.
// returns 0 if list of small primes is exhausted.

void reset(long b);
// resets generator so that the next prime in the sequence
// is the smallest prime >= b.

private:

PrimeSeq(const PrimeSeq&);        // disabled
void operator=(const PrimeSeq&);  // disabled

// auxilliary routines

void start();
void shift(long);

};




/**************************************************************

                      Input/Output

***************************************************************/

istream& operator>>(istream& s, ZZ& x);  
ostream& operator<<(ostream& s, const ZZ& a); 




/******************************************************************

MaxAllocBlock:

The routines managing ZZVec's and vec_ZZ_p's attempt to
allocate space for ZZ's in contiguous blocks of long's from
the memory allocator (malloc/free);  however, they do not
attempt to allocate more than MaxAllocBlock contiguous long's.

Increasing MaxAllocBlock could increase the external fragmentation
of the free pool, and cause the program to run out of memory.
Decreasing MaxAllocBlock could cause the program to run more slowly,
as it will make more call to the storage allocator.

*******************************************************************/


const long MaxAllocBlock = 10000;



/****************************************************************

    Single-precision modular arithmetic

*****************************************************************/


/*
these routines implement single-precision modular arithmetic.
If n is the modulus, all inputs should be in the range 0..n-1.
The number n itself should be in the range 1..2^{NTL_NBITS}-1.
*/



inline long AddMod(long a, long b, long n)
// return (a+b)%n

{
   long res = a + b;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   return res;
#else
   if (res >= n)
      return res - n;
   else
      return res;
#endif
}

inline long SubMod(long a, long b, long n)
// return (a-b)%n

{
   long res = a - b;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   return res;
#else
   if (res < 0)
      return res + n;
   else
      return res;
#endif
}

inline long NegateMod(long a, long n)
{
   return SubMod(0, a, n);
}


#if (defined(NTL_SINGLE_MUL))


#if (!defined(FAST_INT_MUL))


inline long MulMod(long a, long b, long n)
// return (a*b)%n

{
   double ab;
   long q, res;

   ab = ((double) a) * ((double) b);
   q  = (long) (ab/((double) n));  // q could be off by (+/-) 1
   res = (long) (ab - ((double) q)*((double) n));
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}

/*
The following MulMod takes a fourth argument, ninv,
which is assumed to equal 1/((double) n).
It is usually faster than the above.
*/

inline long MulMod(long a, long b, long n, double ninv)
{
   double ab;
   long q, res;

   ab = ((double) a) * ((double) b);
   q  = (long) (ab*ninv);   // q could be off by (+/-) 1
   res = (long) (ab - ((double) q)*((double) n));
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}

/*
Yet another MulMod.
This time, the 4th argument should be ((double) b)/((double) n).
*/

inline long MulMod2(long a, long b, long n, double bninv)
{
   double ab;
   long q, res;

   ab = ((double) a)*((double) b);
   q = (long) (((double) a)*bninv);
   res = (long) (ab - ((double) q)*((double) n));
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}

inline long MulDivRem(long& qq, long a, long b, long n, double bninv)
{
   double ab;
   long q, res;

   ab = ((double) a)*((double) b);
   q = (long) (((double) a)*bninv);
   res = (long) (ab - ((double) q)*((double) n));
   if (res >= n) {
      res -= n;
      q++;
   } else if (res < 0) {
      res += n;
      q--;
   }

   qq = q;
   return res;
}

#else

inline long MulMod(long a, long b, long n)
// return (a*b)%n

{
   double ab, xx;
   long iab, q, res;

   ab = ((double) a) * ((double) b);
   q  = (long) (ab/((double) n));  // q could be off by (+/-) 1

   xx = ab + 4503599627370496.0;
   NTL_FetchLo(iab, xx);

   res = iab - q*n;

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}

/*
The following MulMod takes a fourth argument, ninv,
which is assumed to equal 1/((double) n).
It is usually faster than the above.
*/

inline long MulMod(long a, long b, long n, double ninv)
{
   double ab, xx;
   long iab, q, res;

   ab = ((double) a) * ((double) b);
   q  = (long) (ab*ninv);   // q could be off by (+/-) 1

   xx = ab + 4503599627370496.0;
   NTL_FetchLo(iab, xx);

   res = iab - q*n;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}

/*
Yet another MulMod.
This time, the 4th argument should be ((double) b)/((double) n).
*/

inline long MulMod2(long a, long b, long n, double bninv)
{
   long q, res;

   q = (long) (((double) a)*bninv);
   res = a*b - q*n;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}


inline long MulDivRem(long& qq, long a, long b, long n, double bninv)
{
   long q, res;

   q = (long) (((double) a)*bninv);
   res = a*b - q*n;
   if (res >= n) {
      res -= n;
      q++;
   } else if (res < 0) {
      res += n;
      q--;
   }

   qq = q;
   return res;
}

#endif




#else

inline long MulMod(long a, long b, long n)
{
   long q, res;

   q  = (long) ((((double) a) * ((double) b)) / ((double) n)); 
   res = a*b - q*n;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}

inline long MulMod(long a, long b, long n, double ninv)
{
   long q, res;

   q  = (long) ((((double) a) * ((double) b)) * ninv); 
   res = a*b - q*n;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}


inline long MulMod2(long a, long b, long n, double bninv)
{
   long q, res;

   q  = (long) (((double) a) * bninv);
   res = a*b - q*n;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (NTL_BITS_PER_LONG-1)) & n;
#else
   if (res >= n)
      res -= n;
   else if (res < 0)
      res += n;
#endif
   return res;
}

inline long MulDivRem(long& qq, long a, long b, long n, double bninv)
{
   long q, res;

   q  = (long) (((double) a) * bninv);
   res = a*b - q*n;
   if (res >= n) {
      res -= n;
      q++;
   } else if (res < 0) {
      res += n;
      q--;
   }

   qq = q;
   return res;
}


#endif


long InvMod(long a, long n);
// computes a^{-1} mod n.  Error is raised if undefined.

long PowerMod(long a, long e, long n);
// computes a^e mod n, e >= 0





#endif

