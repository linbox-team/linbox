

/***********************************************************************

   This software is for research and educational purposes only.

************************************************************************/



#ifndef ZZ__H
#define ZZ__H


/********************************************************

   LIP INTERFACE 

   The class ZZ implements signed, arbitrary length integers.

   See IMPLEMENTATION NOTES below.

**********************************************************/


#include "lip.h"
#include "tools.h"


class ZZ {
public:

long *rep; // This is currently public for "emergeny" situations
           // May be private in future versions.

private:

static ZZ _zero;


public:

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


ZZ(long a) { rep = 0; zintoz(a, &rep); }


void operator=(const ZZ& a) 

{ zcopy(a.rep, &rep); }

void operator=(long a) { zintoz(a, &rep); }


~ZZ()

{ zfree(&rep); }

void kill()
// force the space held by this ZZ to be released.
// The value then becomes 0.

{ zfree(&rep); }

void SetSize(long k)
// pre-allocates space for k-digit numbers;  does not change the value.

{ zsetlength(&rep, k); }

long size() const; 
// returns the number of (ZZ_NBIT-bit) digits of |a|; the size of 0 is 0.

static const ZZ& zero() { return _zero; }


inline ZZ(const ZZ& a, const ZZ& b, INIT_ADD_TYPE);
inline ZZ(const ZZ& a, const ZZ& b, INIT_SUB_TYPE);
inline ZZ(const ZZ& a, const ZZ& b, INIT_MUL_TYPE);
inline ZZ(const ZZ& a, const ZZ& b, INIT_DIV_TYPE);
inline ZZ(const ZZ& a, const ZZ& b, INIT_REM_TYPE);
inline ZZ(const ZZ& a, INIT_NEG_TYPE);

inline ZZ(INIT_VAL_TYPE, const char *);

};




inline void clear(ZZ& x)
// x = 0

   { zzero(&x.rep); }

inline void set(ZZ& x)
// x = 1

   { zsetlength(&x.rep, 1); x.rep[0] = 1; x.rep[1] = 1; }


inline void swap(ZZ& x, ZZ& y)
// swap the values of x and y (swaps pointers only)

   { zswap(&x.rep, &y.rep); }



/**********************************************************

   Conversion routines.

***********************************************************/



inline void operator<<(ZZ& x, long a)
//  x = a

{ zintoz(a, &x.rep); }



inline void operator<<(ZZ& x, int a) 
// x = a

{ zintoz((long) a, &x.rep); }

void operator<<(ZZ& x, const char *s);

inline ZZ::ZZ(INIT_VAL_TYPE, const char *s)
{  rep = 0; *this << s; }



inline void operator<<(ZZ& x, double a)
// x = floor(a);

{ zdoubtoz(a, &x.rep); }

inline long Long(const ZZ& a)  // retained for backward compatibility
// return a, no overflow check

   { return ztoint(a.rep); }


/* these are provided for notational consistency */

inline void operator<<(long& x, const ZZ& a)
// x = a, no overflow check

  { x = ztoint(a.rep); }

inline void operator<<(int& x, const ZZ& a)
// x = a, no overflow check

  { x = (int) ztoint(a.rep); }


inline double Double(const ZZ& a)
// return a, no overflow check

   { return zdoub(a.rep); }

/* these are provided for notational consistency */

inline void operator<<(double& x, const ZZ& a)
//  x = a, no overflow check

   { x = zdoub(a.rep); }


inline void operator<<(float& x, const ZZ& a)
//  x = a, no overflow check

   { x = (float) zdoub(a.rep); }

double log(const ZZ& a);


// The following is provided to disable very wierd conversions here
// and also in RR, xdouble, and quad_float.

class disable_wierd_conversions_1 {

public:

disable_wierd_conversions_1() {}

disable_wierd_conversions_1(int) {}
disable_wierd_conversions_1(long) {}
disable_wierd_conversions_1(float) {}
disable_wierd_conversions_1(double) {}

};

long Long(const disable_wierd_conversions_1&);
double Double(const disable_wierd_conversions_1&);


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

inline long compare(const ZZ& a, long b)
   { return zscompare(a.rep, b); }
inline long operator==(const ZZ& a, long b)
   { return zscompare(a.rep, b) == 0; }
inline long operator!=(const ZZ& a, long b)
   { return zscompare(a.rep, b) != 0; }
inline long operator<(const ZZ& a, long b)
   { return zscompare(a.rep, b) < 0; }
inline long operator>(const ZZ& a, long b)
   { return zscompare(a.rep, b) > 0; }
inline long operator<=(const ZZ& a, long b)
   { return zscompare(a.rep, b) <= 0; }
inline long operator>=(const ZZ& a, long b)
   { return zscompare(a.rep, b) >= 0; }

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
// z = a - b;  assumes a >= b >= 0.

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

void sub(ZZ& x, const ZZ& a, long b);




/*******************************************************

                 Multiplication.

********************************************************/

inline void mul(ZZ& x, const ZZ& a, const ZZ& b)
// x = a * b

   { zmul(a.rep, b.rep, &x.rep); }


inline void sqr(ZZ& x, const ZZ& a)
// x = a*a

   { zsq(a.rep, &x.rep); }

/* single-precision versions */

inline void mul(ZZ& x, const ZZ& a, long b)
   { zsmul(a.rep, b, &x.rep); }


void MultiMul(ZZ& x, long n, const ZZ* a, const long* b, long size);
// x = sum_{i=0}^{n-1} a[i]*b[i].
// It is assumed that a is nonnegative, and that
// 0 <= b[i] < 2^{ZZ_NBITS} for each i.
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


/* single-precision versions */

inline long DivRem(ZZ& q, const ZZ& a, long b)
   { return zsdiv(a.rep, b, &q.rep); } 
inline void div(ZZ& q, const ZZ& a, long b)
   { (void) zsdiv(a.rep, b, &q.rep); }
inline long rem(const ZZ& a, long b)
   { return zsmod(a.rep, b); }


long divide(ZZ& q, const ZZ& a, long b);
// if b | a, sets q = a/b and returns 1; otherwise returns 0.


long divide(const ZZ& a, long b);
// if b | a, returns 1; otherwise returns 0.

inline long divide(long a, const ZZ& b) { return divide(ZZ(a), b); }



#if (defined(SINGLE_MUL))

inline void MultiRem(long* r, long n, const ZZ& a, 
                     const long* d, double **tbl)

   { zmultirem2(a.rep, n, (long *) d, tbl, r); }

#elif (defined(TBL_REM))

inline void MultiRem(long* r, long n, const ZZ& a, 
                     const long* d, long **tbl)

   { zmultirem3(a.rep, n, (long *) d, tbl, r); }

#else 

inline void MultiRem(long* r, long n, const ZZ& a, const long* d)

   { zmultirem(a.rep, n, (long *) d, r); }

#endif

// computes r[i] = a % d[i], i = 0..n-1.
// assumes a >= 0 and 0 < d[i] < 2^{ZZ_NBITS}
// This must be fast, as it is critical in chinese remaindering.


/**********************************************************

                        GCD's

***********************************************************/


inline void GCD(ZZ& d, const ZZ& a, const ZZ& b)
// d = gcd(a, b)

   { zgcd(a.rep, b.rep, &d.rep); }


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
// x = (a << k), k < 0 implies right-shift

   { zlshift(a.rep, k, &x.rep); }


inline void RightShift(ZZ& x, const ZZ& a, long k)
// x = (a >> k), k < 0 implies left-shift

   { zrshift(a.rep, k, &x.rep); }

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

inline void LowBits(ZZ& x, const ZZ& a, long k)
// puts k low order bits of |a| into x

   { zlowbits(a.rep, k, &x.rep); }

inline long LowBits(const ZZ& a, long k)
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

inline void and(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| AND |b|

   { zand(a.rep, b.rep, &x.rep); }

inline void or(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| OR |b|

   { zor(a.rep, b.rep, &x.rep); }

inline void xor(ZZ& x, const ZZ& a, const ZZ& b)
// x = |a| XOR |b|

   { zxor(a.rep, b.rep, &x.rep); }


// single-precision versions

long NumBits(long a);

long bit(long a, long k);

long NextPowerOfTwo(long m);
// returns least nonnegative k such that 2^k >= m


/***********************************************************

                  Psuedo-random Numbers

************************************************************/


inline void SetSeed(const ZZ& s)
// initialize random number generator

   { zrstart(s.rep); }

inline void RandomBnd(ZZ& x, const ZZ& n)
// x = "random number" in the range 0..n-1, or 0  if n <= 0

   { zrandomb(n.rep, &x.rep); } 


inline void RandomLen(ZZ& x, long NumBits)
// x = "random number" with precisely NumBits bits.

   { zrandoml(NumBits, &x.rep); }





// single-precision version of the above

inline long RandomBnd(long n)
   { return zrandom(n); }


long RandomLen(long l);


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




/************************************************************

                      Primality Testing 

*************************************************************/


long ProbPrime(const ZZ& n, long NumTrials = 10);
// tests if n is prime;  performs a little trial division,
// followed by a single-precision MillerWitness test, followed by
// up to NumTrials general MillerWitness tests.

long MillerWitness(const ZZ& n, const ZZ& w);
long MillerWitness(const ZZ& n, long w);
// Tests if w is a witness to primality a la Miller.
// Assumption: n is odd and positive, 0 <= w < n.

void RandomPrime(ZZ& n, long l, long NumTrials=10);
// n =  random l-bit prime

void NextPrime(ZZ& n, const ZZ& m, long NumTrials=10);
// n = smallest prime >= m.

// single-precision versions

long ProbPrime(long n, long NumTrials = 10);

long RandomPrime(long l, long NumTrials=10);

long NextPrime(long l, long NumTrials=10);


/************************************************************

                      Exponentiation

*************************************************************/

inline void power(ZZ& x, const ZZ& a, long e)
// x = a^e, e >= 0

   { zexp(a.rep, e, &x.rep); }

inline void power(ZZ& x, long a, long e)
// x = a^e, e >= 0

   { zexps(a, e, &x.rep); } 



void power2(ZZ& x, long e);


/*************************************************************

                       Square Roots

**************************************************************/




inline void SqrRoot(ZZ& x, const ZZ& a)
// x = [a^{1/2}], a >= 0

{
   zsqrt(a.rep, &x.rep);
}


inline long SqrRoot(long a) { return zsqrts(a); }
// single-precision version


/************************************************************

                   Some syntactic sugar

************************************************************/

inline ZZ::ZZ(const ZZ& a, const ZZ& b, INIT_ADD_TYPE)
{
   rep = 0;
   add(*this, a, b);
}

inline ZZ operator+(const ZZ& a, const ZZ& b)
{
   return ZZ(a, b, INIT_ADD);
}

inline ZZ::ZZ(const ZZ& a, const ZZ& b, INIT_SUB_TYPE)
{
   rep = 0;
   sub(*this, a, b);
}

inline ZZ operator-(const ZZ& a, const ZZ& b)
{
   return ZZ(a, b, INIT_SUB);
}

inline ZZ::ZZ(const ZZ& a, const ZZ& b, INIT_MUL_TYPE)
{
   rep = 0;
   mul(*this, a, b);
}

inline ZZ operator*(const ZZ& a, const ZZ& b)
{
   return ZZ(a, b, INIT_MUL);
}


inline ZZ::ZZ(const ZZ& a, const ZZ& b, INIT_DIV_TYPE)
{
   rep = 0;
   div(*this, a, b);
}

inline ZZ operator/(const ZZ& a, const ZZ& b)
{
   return ZZ(a, b, INIT_DIV);
}


inline ZZ::ZZ(const ZZ& a, const ZZ& b, INIT_REM_TYPE)
{
   rep = 0;
   rem(*this, a, b);
}

inline ZZ operator%(const ZZ& a, const ZZ& b)
{
   return ZZ(a, b, INIT_REM);
}


inline void operator +=(ZZ& x, const ZZ& a)
{
   add(x, x, a);
}

inline void operator -=(ZZ& x, const ZZ& a)
{
   sub(x, x, a);
}

inline void operator *=(ZZ& x, const ZZ& a)
{
   mul(x, x, a);
}

inline void operator /=(ZZ& x, const ZZ& a)
{
   div(x, x, a);
}


inline void operator %=(ZZ& x, const ZZ& a)
{
   rem(x, x, a);
}

inline ZZ::ZZ(const ZZ& a, INIT_NEG_TYPE)
{
   rep = 0;
   negate(*this, a);
}

inline ZZ operator-(const ZZ& a)
{
   return ZZ(a, INIT_NEG);
}

inline void operator ++(ZZ& x)
{
   add(x, x, 1);
}

inline void operator ++(ZZ& x, int)
{
   add(x, x, 1);
}

inline void operator --(ZZ& x)
{
   add(x, x, -1);
}

inline void operator --(ZZ& x, int)
{
   add(x, x, -1);
}


ZZ abs(const ZZ& a);



/*************************************************************


                    Small Prime Generation


*************************************************************/


// primes are generated in sequence, starting at 2, 
// and up until (2*PRIME_BND+1)^2, which is less than ZZ_RADIX.

#if (ZZ_NBITS > 30)
#define PRIME_BND ((1L << 14) - 1)
#else
#define PRIME_BND ((1L << (ZZ_NBITS/2-1)) - 1)
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




/***************************************************************

                      Modular Arithmetic

***************************************************************/

// The following routines perform arithmetic mod n, n positive.
// All args (other than exponents) are assumed to be in the range 0..n-1.
// ALIAS RESTRICTION: in all of these routines, it is 
//                    assumed that n is not aliased by 
//                    any of the outputs.




inline void AddMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
// x = (a+b)%n

   { zaddmod(a.rep, b.rep, n.rep, &x.rep); }

inline void SubMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
// x = (a-b)%n

   { zsubmod(a.rep, b.rep, n.rep, &x.rep); }

inline void NegateMod(ZZ& x, const ZZ& a, const ZZ& n)
// x = -a % n

   { zsubmod(0, a.rep, n.rep, &x.rep); }

void AddMod(ZZ& x, const ZZ& a, long b, const ZZ& n);
void AddMod(ZZ& x, long a, const ZZ& b, const ZZ& n);
void SubMod(ZZ& x, const ZZ& a, long b, const ZZ& n);
void SubMod(ZZ& x, long a, const ZZ& b, const ZZ& n);

inline void MulMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
// x = (a*b)%n

   { zmulmod(a.rep, b.rep, n.rep, &x.rep); }

inline void MulMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
// x = (a*b)%n

   { zsmulmod(a.rep, b, n.rep, &x.rep); }

inline void MulMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
   { MulMod(x, b, a, n); }


inline void SqrMod(ZZ& x, const ZZ& a, const ZZ& n)
// x = a^2 % n

   { zsqmod(a.rep, n.rep, &x.rep); }

inline void InvMod(ZZ& x, const ZZ& a, const ZZ& n)
// x = a^{-1} mod n, 0 <= x < n
// error is raised occurs if inverse not defined

   { zinvmod(a.rep, n.rep, &x.rep); }


inline long InvModStatus(ZZ& x, const ZZ& a, const ZZ& n)
// if gcd(a,b) = 1, then ReturnValue = 0, x = a^{-1} mod n
// otherwise, ReturnValue = 1, x = gcd(a, n)

  { return zinv(a.rep, n.rep, &x.rep); }


void PowerMod(ZZ& x, const ZZ& a, const ZZ& e, const ZZ& n);
// x = a^e % n, e >= 0


// single-precision versions

inline void PowerMod(ZZ& x, long a, const ZZ& e, const ZZ& n)
   { zexpmods(a, e.rep, n.rep, &x.rep); }



/*************************************************************

   Jacobi symbol and modular squre roots

**************************************************************/


long Jacobi(const ZZ& a, const ZZ& n);
//  compute Jacobi symbol of a and n;
//  assumes 0 <= a < n, n odd

void SqrRootMod(ZZ& x, const ZZ& a, const ZZ& n);
//  computes square root of a mod n;
//  assumes n is an odd prime, and that a is a square mod n





/**************************************************************

                      Input/Output

***************************************************************/

istream& operator>>(istream& s, ZZ& x);  
ostream& operator<<(ostream& s, const ZZ& a); 




/******************************************************************

MaxAllocBlock:

The routines managing ZZVec's and vector(ZZ_p)'s attempt to
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
The number n itself should be in the range 1..2^{ZZ_NBITS}-1.
*/



inline long AddMod(long a, long b, long n)
// return (a+b)%n

{
   long res = a + b;
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   return res;
#else
   if (res < 0)
      return res + n;
   else
      return res;
#endif
}


#if (defined(SINGLE_MUL))


#if (!defined(FAST_INT_MUL))


inline long MulMod(long a, long b, long n)
// return (a*b)%n

{
   double ab;
   long q, res;

   ab = ((double) a) * ((double) b);
   q  = (long) (ab/((double) n));  // q could be off by (+/-) 1
   res = (long) (ab - ((double) q)*((double) n));
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
   FetchLo(iab, xx);

   res = iab - q*n;

#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
   FetchLo(iab, xx);

   res = iab - q*n;
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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
#if (ZZ_ARITH_RIGHT_SHIFT && defined(AVOID_BRANCHING))
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
   res -= n;
   res += (res >> (ZZ_BITS_PER_LONG-1)) & n;
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

