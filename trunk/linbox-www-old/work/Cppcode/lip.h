
#ifndef lip__H
#define lip__H


/***********************************************************************

   This software is for research and educational purposes only.

************************************************************************/


/***********************************************************************

   This Long Integer Package is derived from that developed
   by Arjen Lenstra.   However, there are many changes,
   and it is generally intended only to be used through the
   C++ interface defined in the ZZ module.

   Outputs may alias inputs, except where explicitly noted.

************************************************************************/


typedef long * verylong;

#include "mach_desc.h"

#define ZZ_NBITSH          (ZZ_NBITS>>1)
#define ZZ_RADIX           (1L<<ZZ_NBITS)
#define ZZ_RADIXM          (ZZ_RADIX-1)
#define ZZ_RADIXROOT       (1L<<ZZ_NBITSH)
#define ZZ_RADIXROOTM      (ZZ_RADIXROOT-1)

#define ZZ_FRADIX ((double) ZZ_RADIX)
#define ZZ_FRADIX_INV  (((double) 1.0)/((double) ZZ_RADIX))



#if (defined(SINGLE_MUL) && !ZZ_SINGLE_MUL_OK)
#undef SINGLE_MUL
#endif

#if (defined(SINGLE_MUL))

#undef ZZ_NBITS
#define ZZ_NBITS (26)

/****************************************************************

The following macros extract the two words of a double,
getting around the type system.
This is only used in the SINGLE_MUL strategy.

*****************************************************************/

#if (ZZ_DOUBLES_LOW_HIGH)
#define ZZ_LO_WD 0
#define ZZ_HI_WD 1
#else
#define ZZ_LO_WD 1
#define ZZ_HI_WD 0
#endif


typedef union { double d; unsigned long rep[2]; } d_or_rep;

#define FetchHiLo(hi,lo,x) \
do { \
   d_or_rep ll_xx; \
   ll_xx.d = (x); \
   hi = ll_xx.rep[ZZ_HI_WD]; \
   lo = ll_xx.rep[ZZ_LO_WD]; \
} while (0)


#define FetchLo(lo,x)  \
do {  \
   d_or_rep ll_xx;  \
   ll_xx.d = x;  \
   lo = ll_xx.rep[ZZ_LO_WD];  \
} while (0) 

#endif


/**********************************************************************

   A multiprecision integer is represented as a pointer to long.
   Let x be such a pointer.
   x = 0 represents 0.
   Otherwise, let n = abs(x[0]) and s = sign(x[0]);
   the integer represented is then:

      s*(x[1] + x[2]*ZZ_RADIX + ... + x[n]*ZZ_RADIX^{n-1}),

   where x[n] != 0, unless n = s = 1.
   Notice that the number zero can be represented in precisely 2 ways,
   either as a null pointer, or as x[0] = 1 and x[1] = 0.

   Storage is generally managed via zsetlength and zfree.
   x[-1] = (k << 1) | b, where k is the maximum value of n allocated,
   and b is a bit that is set is the space is managed by some
   mechanism other than zsetlength and zfree.

************************************************************************/


#if (defined(__cplusplus) && !defined(CPLUSPLUS_ONLY))
extern "C" {
#endif


/***********************************************************************

   Basic Functions

***********************************************************************/
    


    void zsadd(verylong a, long d, verylong *b);
       /* *b = a + d */

    void zadd(verylong a, verylong b, verylong *c);
       /*  *c = a + b */

    void zsub(verylong a, verylong b, verylong *c);
       /* *c = a - b */

    void zsubpos(verylong a, verylong b, verylong *c);
       /* *c = a - b; assumes a >= b >= 0 */

    void zsmul(verylong a, long d, verylong *b);
       /* *b = d * a */

    void zmul(verylong a, verylong b, verylong *c);
       /* *c = a * b */

    void zsq(verylong a, verylong *c);
       /* *c = a * a */

    long zsdiv(verylong a, long b, verylong *q);
       /* (*q) = floor(a/b) and a - floor(a/b)*(*q) is returned;
          error is raised if b == 0;
          if b does not divide a, then sign(*q) == sign(b) */

    void zdiv(verylong a, verylong b, verylong *q, verylong *r);
       /* (*q) = floor(a/b) and (*r) = a - floor(a/b)*(*q);
          error is raised if b == 0;
          if b does not divide a, then sign(*q) == sign(b) */

    void zmultirem(verylong a, long n, long* dd, long* rr);
    void zmultirem2(verylong a, long n, long* dd, double **tbl, long* rr);
       /* rr[i] = a % dd[i], i = 0..n-1;
          assumes a >= 0, 0 < dd[i] < ZZ_RADIX
          zmultirem2 takes an extra argument, tbl, which contains
          pre-computed residues of powers of RADIX */
    void zmultirem3(verylong a, long n, long* dd, long **tbl, long* rr);
       /* same as above, but tbl has different type */

    long zsfastrem(verylong a, long d);
       /* return a % d;
          assumes a >= 0, 0 < d < ZZ_RADIX */
        

    void zmod(verylong a, verylong b, verylong *r);
       /* same as zdiv, but only remainder is computed */

    long zsmod(verylong a, long d);
       /* same as zsdiv, but only remainder is computed */

    void zquickmod(verylong *r, verylong b);
       /* *r = *r % b;
	  assumes b > 0 and *r >= 0;
	  The division is performed in place (but may sometimes
          cause *r to grow by one digit) */

/********************************************************************

   Shifting and bit manipulation

*********************************************************************/

    void z2mul(verylong n, verylong *a);
       /* *a = 2 * n */

    long z2div(verylong n, verylong *a);
       /* *a = sign(n) * (|n|/2) */ 

    void zlshift(verylong n, long k, verylong *a);
       /* *a = sign(n) * (|n| << k);
          shift is in reverse direction for negative k */

    void zrshift(verylong n, long k, verylong *a);
       /* *a = sign(n) * (|n| >> k);
          shift is in reverse direction for negative k */
    
    long zmakeodd(verylong *n);
       /*
          if (n != 0)
              *n = m;
              return (k such that n == 2 ^ k * m with m odd);
          else
              return (0); 
        */

    long zodd(verylong a);
       /* returns 1 if n is odd and 0 if it is even */

    long zbit(verylong a, long p);
       /* returns p-th bit of a, where the low order bit is indexed by 0;
          p out of range returns 0 */

    long zsetbit(verylong *a, long p);
       /* returns original value of p-th bit of |a|, and replaces
          p-th bit of a by 1 if it was zero;
          error if p < 0 */

    long zswitchbit(verylong *a, long p);
       /* returns original value of p-th bit of |a|, and switches
          the value of p-th bit of a;
          p starts counting at 0;
          error if p < 0 */


     void zlowbits(verylong a, long k, verylong *b);
        /* places k low order bits of |a| in b */ 

     long zslowbits(verylong a, long k);
        /* returns k low order bits of |a| */

    long zweights(long a);
        /* returns Hamming weight of |a| */

    long zweight(verylong a);
        /* returns Hamming weight of |a| */

    void zand(verylong a, verylong b, verylong *c);
        /* c gets bit pattern `bits of |a|` and `bits of |b|` */

    void zor(verylong a, verylong b, verylong *c);
        /* c gets bit pattern `bits of |a|` inclusive or `bits of |b|` */

    void zxor(verylong a, verylong b, verylong *c);
        /* c gets bit pattern `bits of |a|` exclusive or `bits of |b|` */




/************************************************************************

   Comparison

*************************************************************************/

    long zcompare(verylong a, verylong b);
       /*
          if (a > b)
              return (1);
          if (a == b)
              return (0);
          if (a < b)
              return (-1);
         */

    long zscompare(verylong a, long b);
       /* single-precision version of the above */

    long ziszero (verylong a);
       /* test for 0 */


    long zsign(verylong a);
       /* 
          if (a > 0)
              return (1);
          if (a == 0)
              return (0);
          if (a < 0)
              return (-1);
        */

    void zabs(verylong *a);
       /* *a = |a| */

    void znegate(verylong *a);
       /* *a = -a */

    void zcopy(verylong a, verylong *b);
       /* *b = a */

    void zswap(verylong *a, verylong *b);
       /* swap a and b (by swaping pointers) */

    long z2log(verylong a);
       /* number of bits in |a|; returns 0 if a = 0 */

    long z2logs(long a);
        /* single-precision version of the above */


/********************************************************************

   Conversion

*********************************************************************/
        
    void zzero(verylong *a);
       /* *a = 0;  space is allocated */

    void zone(verylong *a);
       /* *a = 1 */

    void zintoz(long d, verylong *a);
       /* *a = d */

    long ztoint(verylong a);
       /* converts a to a long;  overflow results in value
          mod 2^{ZZ_BITS_PER_LONG}. */


    double zdoub(verylong n);
       /* converts a to a double;  no overflow check */

    void zdoubtoz(double a, verylong *x);
       /* x = floor(a);  
          assumes |floor(a)| < 2^{ZZ_DOUBLE_PRECISION-1} */
    



/************************************************************************

   Square roots

*************************************************************************/


    long zsqrts(long n);
       /* return floor(sqrt(n));  error raised in n < 0 */

    void zsqrt(verylong n, verylong *r);
       /* *r =  floor(sqrt(n));  error raised in n < 0 */

/*********************************************************************
 
    Exponentiation
 
**********************************************************************/

   void zexp(verylong a, long e, verylong *b);
       /* *b = a^e;  error raised if e < 0 */

   void zexps(long a, long e, verylong *b);
       /* *b = a^e;  error raised if e < 0 */
       

/*********************************************************************

   Modular Arithmetic

   Addition, subtraction, multiplication, squaring division, inversion,
   and exponentiation modulo a positive modulus n, where all operands
   (except for the exponent in exponentiation) and results are in the
   range [0, n-1].   

   ALIAS RESTRICTION:  output parameters should not alias n

***********************************************************************/

    void zaddmod(verylong a, verylong b, verylong n, verylong *c);
       /* *c = (a + b) % n */

    void zsubmod(verylong a, verylong b, verylong n, verylong *c);
       /* *c = (a - b) % n */

    void zsmulmod(verylong a, long b, verylong n, verylong *c);
       /* *c = (a * b) % n */

    void zmulmod(verylong a, verylong b, verylong n, verylong *c);
       /* *c = (a * b) % n */

    void zsqmod(verylong a, verylong n, verylong *c);
       /* *c = (a ^ 2) % n */

    void zinvmod(verylong a, verylong n, verylong *c);
       /* *c = (1 / a) % n; error raised if gcd(b, n) != 1 */


    void zsexpmod(verylong a, long e, verylong n, verylong *b);
    void zsexpmods(long a, long e, verylong n, verylong *b);
    void zexpmod(verylong a, verylong e, verylong n, verylong *b);
    void zexpmods(long a, verylong e, verylong n, verylong *b);
       /* *b = (a ^ e) % n; error raised if e < 0 */



/**************************************************************************

   Euclidean Algorithms

***************************************************************************/
    void zgcd(verylong m1, verylong m2, verylong *r);
       /* *r = greatest common divisor of m1 and m2; 
          uses binary gcd algorithm */


    void zexteucl(verylong a, verylong *xa,
                 verylong b, verylong *xb,
                 verylong *d);
       /*
          *d = a * *xa + b * *xb = gcd(a, b);
          sets *d, *xa and *xb given a and b;
          uses Lehmer`s trick
        */


    long zinv(verylong a, verylong b, verylong *c);
       /*
          if (a and b coprime)
          {
              *c = inv; 
              return(0);
          }
          else
          {
              *c = gcd(a, b);
              return(1);
          }
          
          where inv is such that (inv * a)  == 1 mod b;
          error raised if a < 0 or b <= 0
        */

/**************************************************************************

   Random number generation

   These algorithms produce pseudo-random numbers that
   can be used for probabilistic algorithms, but certainly
   should not be used for cryptographic purposes

***************************************************************************/


    void zrstart(verylong s);
       /* initialize the random number generator with seed s;
          if this is not called, the default initial value is 7157891. */

    long zrandom(long b);
       /* return value is a pseudo-random integer in [0, b) */

    void zrandomb(verylong bnd, verylong *a);
       /*
          sets a to a random integer in the interval [0, bnd),
          a gets zero if interval is empty.
         
          zrandomb uses zrandom to generate large enough random number,
          which is modded out by bnd to get the resulting a
        */

    void zrandoml(long length, verylong *a);
       /* sets a to a random integer of precisely length bits */

        
/**********************************************************************

    Storage Allocation

    These routines use malloc and free.

***********************************************************************/


    void zsetlength(verylong *v, long len);
       /* Allocates enough space to hold a len-digit number,
          where each digit has ZZ_NBITS bits.
          If space must be allocated, space for one extra digit
          is always allocated. */

    void zfree(verylong *x);
       /* Free's space held by x, and sets x back to 0. */


/*******************************************************************

    Special routines

********************************************************************/



    void zaddmul(long ams, long *ama, long *amb);
       /* internal use in module ZZ only */



#if (defined(__cplusplus) && !defined(CPLUSPLUS_ONLY))
}
#endif



#endif
