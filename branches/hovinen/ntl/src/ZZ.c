

#include <NTL/ZZ.h>
#include <NTL/tools.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include <NTL/vec_ZZ.h>


const ZZ& ZZ::zero()
{
   static ZZ z;
   return z;
}


const ZZ& ZZ_expo(long e)
{
   ZZ expo_helper;
   conv(expo_helper, e);
   return expo_helper;
}



long ZZ::size() const
{
   if (!rep || (rep[0] == 1 && rep[1] == 0))
      return 0;
   else if (rep[0] < 0)
      return -rep[0];
   else 
      return rep[0];
}

long digit(const ZZ& a, long i)
{
   verylong rep = a.rep;

   if (i < 0 || !rep) return 0;

   long sa = rep[0];
   if (sa < 0) sa = -sa;
   if (i >= sa) return 0;
   return rep[i+1];
}


long IsOne(const ZZ& a)
{
   return a.rep != 0 && a.rep[0] == 1 && a.rep[1] == 1;
}


void AddMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
{
   ZZ B;
   conv(B, b);
   AddMod(x, a, B, n);
}


void SubMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
{
   ZZ B;
   conv(B, b);
   SubMod(x, a, B, n);
}

void SubMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
{
   ZZ A;
   conv(A, a);
   SubMod(x, A, b, n);
}


static
void SPPowerMod(ZZ& x, long a, const ZZ& e, const ZZ& n)
{
   if (IsZero(e)) {
      set(x);
      return;
   }

   ZZ res(INIT_SIZE, n.size());

   long k = NumBits(e);
   long i;

   set(res);

   for (i = k - 1; i >= 0; i--) {
      SqrMod(res, res, n);
      if (bit(e, i))
         MulMod(res, res, a, n);
   }

   if (e < 0) InvMod(res, res, n);

   x = res;
}

static
void TwoPowerMod(ZZ& x, const ZZ& e, const ZZ& n)
{
   if (IsZero(e)) {
      set(x);
      return;
   }

   ZZ res(INIT_SIZE, n.size());

   long k = NumBits(e);
   long i;

   set(res);

   for (i = k - 1; i >= 0; i--) {
      SqrMod(res, res, n);
      if (bit(e, i))
         AddMod(res, res, res, n);
   }

   if (e < 0) InvMod(res, res, n);

   x = res;
}


// ****** input and output

static long iodigits = 0;
static long ioradix = 0;

// iodigits is the greatest integer such that 10^{iodigits} < NTL_RADIX
// ioradix = 10^{iodigits}

static void InitZZIO()
{
   long x;

   x = (NTL_RADIX-1)/10;
   iodigits = 0;
   ioradix = 1;

   while (x) {
      x = x / 10;
      iodigits++;
      ioradix = ioradix * 10;
   }

   if (iodigits <= 0) Error("problem with I/O");
}

istream& operator>>(istream& s, ZZ& x)
{
   long c;
   long sign;
   long ndigits;
   long acc;
   ZZ a;

   if (!s) Error("bad ZZ input");

   if (!iodigits) InitZZIO();

   a = 0;

   c = s.peek();
   while (c == ' ' || c == '\n' || c == '\t') {
      s.get();
      c = s.peek();
   }

   if (c == '-') {
      sign = -1;
      s.get();
      c = s.peek();
   }
   else
      sign = 1;

   if (c < '0' || c > '9') Error("bad ZZ input");

   ndigits = 0;
   acc = 0;
   while (c >= '0' && c <= '9') {
      acc = acc*10 + c - '0';
      ndigits++;

      if (ndigits == iodigits) {
         mul(a, a, ioradix);
         add(a, a, acc);
         ndigits = 0;
         acc = 0;
      }

      s.get();
      c = s.peek();
   }

   if (ndigits != 0) {
      long mpy = 1;
      while (ndigits > 0) {
         mpy = mpy * 10;
         ndigits--;
      }

      mul(a, a, mpy);
      add(a, a, acc);
   }

   if (sign == -1)
      negate(a, a);

   x = a;

   return s;
}

struct lstack {
   long top;
   long alloc;
   long *elts;

   lstack() { top = -1; alloc = 0; elts = 0; }
   ~lstack() { }

   long pop() { return elts[top--]; }
   long empty() { return (top == -1); }
   void push(long x);
};

void lstack::push(long x)
{
   if (alloc == 0) {
      alloc = 100;
      elts = (long *) malloc(alloc * sizeof(long));
   }

   top++;

   if (top + 1 > alloc) {
      alloc = 2*alloc;
      elts = (long *) realloc(elts, alloc * sizeof(long));
   }

   if (!elts) {
      Error("out of space in ZZ output");
   }

   elts[top] = x;
}


static
void PrintDigits(ostream& s, long d, long justify)
{
   char *buf = 0;

   if (!buf) {
      buf = (char *) malloc(iodigits);
      if (!buf) Error("out of memory");
   }

   long i = 0;

   while (d) {
      buf[i] = (d % 10) + '0';
      d = d / 10;
      i++;
   }

   if (justify) {
      long j = iodigits - i;
      while (j > 0) {
         s << "0";
         j--;
      }
   }

   while (i > 0) {
      i--;
      s << buf[i];
   }
}
      

   

ostream& operator<<(ostream& s, const ZZ& a)
{
   ZZ b;
   lstack S;
   long r;
   long k;

   if (!iodigits) InitZZIO();

   b = a;

   k = sign(b);

   if (k == 0) {
      s << "0";
      return s;
   }

   if (k < 0) {
      s << "-";
      negate(b, b);
   }

   do {
      r = DivRem(b, b, ioradix);
      S.push(r);
   } while (!IsZero(b));

   r = S.pop();
   PrintDigits(s, r, 0);

   while (!S.empty()) {
      r = S.pop();
      PrintDigits(s, r, 1);
   }
      
   return s;
}


// ******  MultiMul 


void MultiMul(ZZ& x, long n, const ZZ* a, const long* b, long size)
{
   verylong xx, yy;
   long i, sx;

   sx = size+1;

   zsetlength(&x.rep, sx);
   xx = x.rep;

   for (i = 1; i <= sx; i++)
      xx[i] = 0;

   xx++;

   for (i = 0; i < n; i++) {
      yy = a[i].rep;

      if (!yy || !b[i]) continue;

      zaddmul(b[i], xx, yy);
      yy = xx + yy[0];
   
      if ((*yy) >= NTL_RADIX) {
         (*yy) -= NTL_RADIX;
         yy++;
         while ((*yy) == NTL_RADIX-1) {
            *yy = 0;
            yy++;
         }
         (*yy)++;
      }
   }

   xx--;
   while (sx > 1 && xx[sx] == 0) sx--;
   xx[0] = sx;
}


long GCD(long a, long b)
{
   long u, v, t, x;

   if (a < 0)
      a = -a;

   if (b < 0)
      b = -b;

   if (a < 0 || b < 0)
      Error("GCD: integer overflow");

   if (b==0)
      x = a;
   else {
      u = a;
      v = b;
      do {
         t = u % v;
         u = v; 
         v = t;
      } while (v != 0);

      x = u;
   }

   return x;
}

         

void XGCD(long& d, long& s, long& t, long a, long b)
{
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      b = -b;
      bneg = 1;
   }

   if (a < 0 || b < 0)
      Error("XGCD: integer overflow");

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d = u;
   s = u1;
   t = v1;
}
   

long InvMod(long a, long n)
{
   long d, s, t;

   XGCD(d, s, t, a, n);
   if (d != 1) Error("InvMod: inverse undefined");
   if (s < 0)
      return s + n;
   else
      return s;
}


long PowerMod(long a, long ee, long n)
{
   long x, y;

   unsigned long e;

   if (ee < 0)
      e = -ee;
   else
      e = ee;

   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = MulMod(x, y, n);
      y = MulMod(y, y, n);
      e = e >> 1;
   }

   if (ee < 0) x = InvMod(x, n);

   return x;
}

long ProbPrime(long n, long NumTests)
{
   long m, x, y, z;
   long i, j, k;

   if (n <= 1) return 0;


   if (n == 2) return 1;
   if (n % 2 == 0) return 0;

   if (n == 3) return 1;
   if (n % 3 == 0) return 0;

   if (n == 5) return 1;
   if (n % 5 == 0) return 0;

   if (n == 7) return 1;
   if (n % 7 == 0) return 0;

   if (n >= NTL_RADIX) {
      return ProbPrime(to_ZZ(n), NumTests);
   }

   m = n - 1;
   k = 0;
   while((m & 1) == 0) {
      m = m >> 1;
      k++;
   }

   // n - 1 == 2^k * m, m odd

   for (i = 0; i < NumTests; i++) {
      x = RandomBnd(n);


      if (x == 0) continue;
      z = PowerMod(x, m, n);
      if (z == 1) continue;
   
      j = 0;
      do {
         y = z;
         z = MulMod(y, y, n);
         j++;
      } while (!(j == k || z == 1));

      if (z != 1 || y !=  n-1) return 0;
   }

   return 1;
}


long MillerWitness(const ZZ& n, const ZZ& x)
{
   ZZ m, y, z;
   long j, k;

   if (x == 0) return 0;

   add(m, n, -1);
   k = MakeOdd(m);
   // n - 1 == 2^k * m, m odd

   PowerMod(z, x, m, n);
   if (z == 1) return 0;

   j = 0;
   do {
      y = z;
      SqrMod(z, y, n);
      j++;
   } while (!(j == k || z == 1));

   if (z != 1) return 1;
   add(y, y, 1);
   if (y != n) return 1;
   return 0;
}


long ProbPrime(const ZZ& n, long NumTrials)
{
   if (n <= 1) return 0;

   if (n.size() == 1) {
      return ProbPrime(to_long(n), NumTrials);
   }

   long bn = NumBits(n);
   long wn = (bn+NTL_NBITS-1)/NTL_NBITS;
   long fn = wn/4+1;

   long prime_bnd;

   if (NumBits(bn) + NumBits(fn) > NTL_NBITS)
      prime_bnd = NTL_RADIX;
   else
      prime_bnd = bn*fn;


   PrimeSeq s;
   long p;

   p = s.next();
   while (p && p < prime_bnd) {
      if (rem(n, p) == 0)
         return 0;

      p = s.next();
   }

   ZZ W;
   W = 2;

   // first try W == 2....the exponentiation
   // algorithm runs slightly faster in this case

   if (MillerWitness(n, W))
      return 0;

   long i;

   for (i = 0; i < NumTrials; i++) {
      RandomBnd(W, n);
      if (MillerWitness(n, W)) 
         return 0;
   }

   return 1;
}


void RandomPrime(ZZ& n, long l, long NumTrials)
{
   if (l <= 1)
      Error("RandomPrime: l out of range");

   if (l == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }

   do {
      RandomLen(n, l);
      if (!IsOdd(n)) add(n, n, 1);
   } while (!ProbPrime(n, NumTrials));
}

void NextPrime(ZZ& n, const ZZ& m, long NumTrials)
{
   ZZ x;

   if (m <= 2) {
      n = 2;
      return;
   }

   x = m;

   while (!ProbPrime(x, NumTrials))
      add(x, x, 1);

   n = x;
}

long NextPrime(long m, long NumTrials)
{
   long x;

   if (m <= 2) 
      return 2;

   x = m;

   while (x < NTL_RADIX && !ProbPrime(x, NumTrials))
      x++;

   if (x >= NTL_RADIX)
      Error("NextPrime: no more primes");

   return x;
}



long NextPowerOfTwo(long m)
{
   long k; 
   long n;
   n = 1;
   k = 0;
   while (n < m && n >= 0) {
      n = n << 1;
      k++;
   }

   if (n < 0) Error("NextPowerOfTwo: overflow");

   return k;
}



long NumBits(long a)
{
   unsigned long aa;
   if (a < 0) 
      aa = -a;
   else
      aa = a;

   long k = 0;
   while (aa) {
      k++;
      aa = aa >> 1;
   }

   return k;
}


long bit(long a, long k)
{
   unsigned long aa;
   if (a < 0)
      aa = -a;
   else
      aa = a;

   if (k < 0 || k >= NTL_BITS_PER_LONG) 
      return 0;
   else
      return long((aa >> k) & 1);
}



long divide(ZZ& q, const ZZ& a, const ZZ& b)
{
   ZZ qq, r;

   if (IsZero(b)) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }


   if (IsOne(b)) {
      q = a;
      return 1;
   }

   DivRem(qq, r, a, b);
   if (!IsZero(r)) return 0;
   q = qq;
   return 1;
}

long divide(const ZZ& a, const ZZ& b)
{
   ZZ r;

   if (IsZero(b)) return IsZero(a);
   if (IsOne(b)) return 1;

   rem(r, a, b);
   return IsZero(r);
}

long divide(ZZ& q, const ZZ& a, long b)
{
   ZZ qq;

   if (!b) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }

   if (b == 1) {
      q = a;
      return 1;
   }

   long r = DivRem(qq, a, b);
   if (r) return 0;
   q = qq;
   return 1;
}

long divide(const ZZ& a, long b)
{
   if (!b) return IsZero(a);
   if (b == 1) {
      return 1;
   }

   long r = rem(a,  b);
   return (r == 0);
}
   


long RandomPrime_long(long l, long NumTrials)
{
   if (l <= 1 || l >= NTL_BITS_PER_LONG)
      Error("RandomPrime: l out of range");

   long n;
   do {
      n = RandomLen_long(l);
   } while (!ProbPrime(n, NumTrials));

   return n;
}


PrimeSeq::PrimeSeq()
{
   movesieve = 0;
   movesieve_mem = 0;
   pshift = -1;
   pindex = -1;
   exhausted = 0;
}

PrimeSeq::~PrimeSeq()
{
   if (movesieve_mem)
      free(movesieve_mem);
}

long PrimeSeq::next()
{
   if (exhausted) {
      return 0;
   }

   if (pshift < 0) {
      shift(0);
      return 2;
   }

   for (;;) {
      char *p = movesieve;
      long i = pindex;

      while ((++i) < NTL_PRIME_BND) {
         if (p[i]) {
            pindex = i;
            return pshift + 2 * i + 3;
         }
      }

      long newshift = pshift + 2*NTL_PRIME_BND;

      if (newshift > 2 * NTL_PRIME_BND * (2 * NTL_PRIME_BND + 1)) {
         /* end of the road */
         exhausted = 1;
         return 0;
      }

      shift(newshift);
   }
}

static char *lowsieve = 0;

void PrimeSeq::shift(long newshift)
{
   long i;
   long j;
   long jstep;
   long jstart;
   long ibound;
   char *p;

   if (!lowsieve)
      start();

   pindex = -1;
   exhausted = 0;

   if (newshift < 0) {
      pshift = -1;
      return;
   }

   if (newshift == pshift) return;

   pshift = newshift;

   if (pshift == 0) {
      movesieve = lowsieve;
   } 
   else {
      if (!movesieve_mem) {
         movesieve_mem = (char *) malloc(NTL_PRIME_BND);
         if (!movesieve_mem) 
            Error("out of memory in PrimeSeq");
      }

      p = movesieve = movesieve_mem;
      for (i = 0; i < NTL_PRIME_BND; i++)
         p[i] = 1;

      jstep = 3;
      ibound = pshift + 2 * NTL_PRIME_BND + 1;
      for (i = 0; jstep * jstep <= ibound; i++) {
         if (lowsieve[i]) {
            if (!((jstart = (pshift + 2) / jstep + 1) & 1))
               jstart++;
            if (jstart <= jstep)
               jstart = jstep;
            jstart = (jstart * jstep - pshift - 3) / 2;
            for (j = jstart; j < NTL_PRIME_BND; j += jstep)
               p[j] = 0;
         }
         jstep += 2;
      }
   }
}


void PrimeSeq::start()
{
   long i;
   long j;
   long jstep;
   long jstart;
   long ibnd;
   char *p;

   p = lowsieve = (char *) malloc(NTL_PRIME_BND);
   if (!p)
      Error("out of memory in PrimeSeq");

   for (i = 0; i < NTL_PRIME_BND; i++)
      p[i] = 1;
      
   jstep = 1;
   jstart = -1;
   ibnd = (SqrRoot(2 * NTL_PRIME_BND + 1) - 3) / 2;
   for (i = 0; i <= ibnd; i++) {
      jstart += 2 * ((jstep += 2) - 1);
      if (p[i])
         for (j = jstart; j < NTL_PRIME_BND; j += jstep)
            p[j] = 0;
   }
}

void PrimeSeq::reset(long b)
{
   if (b > (2*NTL_PRIME_BND+1)*(2*NTL_PRIME_BND+1)) {
      exhausted = 1;
      return;
   }

   if (b <= 2) {
      shift(-1);
      return;
   }

   if ((b & 1) == 0) b++;

   shift(((b-3) / (2*NTL_PRIME_BND))* (2*NTL_PRIME_BND));
   pindex = (b - pshift - 3)/2 - 1;
}
 
long Jacobi(const ZZ& aa, const ZZ& nn)
{
   ZZ a, n;
   long t, k;
   long d;

   a = aa;
   n = nn;
   t = 1;

   while (a != 0) {
      k = MakeOdd(a);
      d = trunc_long(n, 3);
      if ((k & 1) && (d == 3 || d == 5)) t = -t;

      if (trunc_long(a, 2) == 3 && (d & 3) == 3) t = -t;
      swap(a, n);
      QuickRem(a, n);
   }

   if (n == 1)
      return t;
   else
      return 0;
}


void SqrRootMod(ZZ& x, const ZZ& aa, const ZZ& nn)
{
   if (aa == 0) {
      x = 0;
      return;
   }

   long i, k;
   ZZ ma, n, t, u, v, e;
   ZZ t1, t2, t3;

   n = nn;
   NegateMod(ma, aa, n);

   // find t such that t^2 - 4*a is not a square

   MulMod(t1, ma, 4, n);
   do {
      RandomBnd(t, n);
      SqrMod(t2, t, n);
      AddMod(t2, t2, t1, n);
   } while (Jacobi(t2, n) != -1);

   // compute u*X + v = X^{(n+1)/2} mod f, where f = X^2 - t*X + a

   add(e, n, 1);
   RightShift(e, e, 1);

   u = 0;
   v = 1;

   k = NumBits(e);

   for (i = k - 1; i >= 0; i--) {
      SqrMod(t1, u, n);
      SqrMod(t2, v, n);
      MulMod(t3, u, v, n);
      MulMod(t3, t3, 2, n);
      MulMod(u, t1, t, n);
      AddMod(u, u, t3, n);
      MulMod(v, t1, ma, n);
      AddMod(v, v, t2, n);

      if (bit(e, i)) {
         MulMod(t1, u, t, n);
         AddMod(t1, t1, v, n);
         MulMod(v, u, ma, n);
         u = t1;
      }

   }

   x = v;
}


// test if a > 0 and -a/2 < g <= a/2
// this is "hand crafted" so as not too waste too much time
// in the CRT routines.

long CRTInRange(const ZZ& gg, const ZZ& aa)
{
   long *g = gg.rep;
   long *a = aa.rep;

   if (!a || a[0] < 0 || (a[0] == 1 && a[1] == 0)) return 0;
   long sa = a[0];

   if (!g) return 1;

   long sg = g[0];
   if (sg == 1 && g[1] == 0) return 1;

   if (sg < 0) sg = -sg;

   if (sa-sg > 1) return 1;

   if (sa-sg < 0) return 0;

   long carry=0;

   if (sa-sg == 1) {
      if (a[sa] > 1) return 1;
      carry = 1;
   }

   long i = sg;
   long diff = 0;
   while (i > 0 && diff == 0) {
      diff = (carry << (NTL_NBITS-1)) + (a[i] >> 1) - g[i];
      carry = (a[i] & 1);
      i--;
   }

   if (diff == 0) {
      if (carry) return 1;
      return (g[0] > 0);
   }
   else
      return (diff > 0);
}


// Chinese Remaindering.
//
// This version in new to v3.7, and is significantly
// simpler and faster than the previous version.
//
// This function takes as input g, a, G, p,
// such that a > 0, 0 <= G < p, and gcd(a, p) = 1.
// It computes a' = a*p and g' such that 
//   * g' = g (mod a);
//   * g' = G (mod p);
//   * -a'/2 < g' <= a'/2.
// It then sets g := g' and a := a', and returns 1 iff g has changed.
//
// Under normal use, the input value g satisfies -a/2 < g <= a/2;
// however, this was not documented or enforced in earlier versions,
// so to maintain backward compatability, no restrictions are placed
// on g.  This routine runs faster, though, if -a/2 < g <= a/2,
// and the first thing the routine does is to make this condition
// hold.
//
// Also, under normal use, both a and p are odd;  however, the routine
// will still work even if this is not so.
//
// The routine is based on the following simple fact.
//
// Let -a/2 < g <= a/2, and let h satisfy
//   * g + a h = G (mod p);
//   * -p/2 < h <= p/2.
// Further, if p = 2*h and g > 0, set
//   g' := g - a h;
// otherwise, set
//   g' := g + a h.
// Then g' so defined satisfies the above requirements.
//
// It is trivial to see that g's satisfies the congruence conditions.
// The only thing is to check that the "balancing" condition
// -a'/2 < g' <= a'/2 also holds.


long CRT(ZZ& gg, ZZ& a, long G, long p)
{
   if (p >= NTL_RADIX) {
      ZZ GG, pp;
      conv(GG, G);
      conv(pp, p);
      return CRT(gg, a, GG, pp);
   }

   long modified = 0;

   ZZ g;

   if (!CRTInRange(gg, a)) {
      modified = 1;
      ZZ a1;
      rem(g, gg, a);
      RightShift(a1, a, 1);
      if (g > a1) sub(g, g, a);
   }
   else
      g = gg;


   long p1;
   p1 = p >> 1;

   long a_inv;
   a_inv = rem(a, p);
   a_inv = InvMod(a_inv, p);

   long h;
   h = rem(g, p);
   h = SubMod(G, h, p);
   h = MulMod(h, a_inv, p);
   if (h > p1)
      h = h - p;

   if (h != 0) {
      modified = 1;
      ZZ ah;
      mul(ah, a, h);

      if (!(p & 1) && g > 0 && (h == p1))
         sub(g, g, ah);
      else
         add(g, g, ah);
   }

   mul(a, a, p);
   gg = g;

   return modified;
}

long CRT(ZZ& gg, ZZ& a, const ZZ& G, const ZZ& p)
{
   long modified = 0;

   ZZ g;

   if (!CRTInRange(gg, a)) {
      modified = 1;
      ZZ a1;
      rem(g, gg, a);
      RightShift(a1, a, 1);
      if (g > a1) sub(g, g, a);
   }
   else
      g = gg;


   ZZ p1;
   RightShift(p1, p, 1);

   ZZ a_inv;
   rem(a_inv, a, p);
   InvMod(a_inv, a_inv, p);

   ZZ h;
   rem(h, g, p);
   SubMod(h, G, h, p);
   MulMod(h, h, a_inv, p);
   if (h > p1)
      sub(h, h, p);

   if (h != 0) {
      modified = 1;
      ZZ ah;
      mul(ah, a, h);

      if (!IsOdd(p) && g > 0 &&  (h == p1))
         sub(g, g, ah);
      else
         add(g, g, ah);
   }

   mul(a, a, p);
   gg = g;

   return modified;
}



void sub(ZZ& x, const ZZ& a, long b)
{
   ZZ B;
   conv(B, b);
   sub(x, a, B);
}

void sub(ZZ& x, long a, const ZZ& b)
{
   ZZ A;
   conv(A, a);
   sub(x, A, b);
}


void power2(ZZ& x, long e)
{
   if (e < 0) Error("power2: negative exponent");
   set(x);
   LeftShift(x, x, e);
}


double log(const ZZ& a)
{
   double log_radix = log(NTL_RADIX);

   if (sign(a) <= 0)
      Error("log(ZZ): argument <= 0");

   long *p = a.rep;

   long i = p[0];

   long prec = 1;
   double x = p[i];
   i--;

   while (i > 0 && prec <= NTL_DOUBLE_PRECISION) {
      x = x * NTL_RADIX + double(p[i]);
      i--;
      prec += NTL_NBITS;
   }


   return log(x) + i*log_radix;
}

   
   
void conv(ZZ& x, const char *s)
{
   long c;
   long sign;
   long ndigits;
   long acc;
   long i = 0;

   ZZ a;

   if (!s) Error("bad ZZ input");

   if (!iodigits) InitZZIO();

   a = 0;

   c = s[i];
   while (c == ' ' || c == '\n' || c == '\t') {
      i++;
      c = s[i];
   }

   if (c == '-') {
      sign = -1;
      i++;
      c = s[i];
   }
   else
      sign = 1;

   if (c < '0' || c > '9') Error("bad ZZ input");

   ndigits = 0;
   acc = 0;
   while (c >= '0' && c <= '9') {
      acc = acc*10 + c - '0';
      ndigits++;

      if (ndigits == iodigits) {
         mul(a, a, ioradix);
         add(a, a, acc);
         ndigits = 0;
         acc = 0;
      }

      i++;
      c = s[i];
   }

   if (ndigits != 0) {
      long mpy = 1;
      while (ndigits > 0) {
         mpy = mpy * 10;
         ndigits--;
      }

      mul(a, a, mpy);
      add(a, a, acc);
   }

   if (sign == -1)
      negate(a, a);

   x = a;
}


static
long OptWinSize(long n)
// finds k that minimizes n/(k+1) + 2^{k-1}

{
   long k;
   double v, v_new;


   v = n/2.0 + 1.0;
   k = 1;

   for (;;) {
      v_new = n/(double(k+2)) + double(1L << k);
      if (v_new >= v) break;
      v = v_new;
      k++;
   }

   return k;
}
      


void PowerMod(ZZ& h, const ZZ& g, const ZZ& e, const ZZ& F)
// h = g^e mod f using "sliding window" algorithm

// remark: the notation (h, g, e, F) is strange, because I
// copied the code from BB.c.

{

   if (sign(g) < 0 || g >= F || F <= 1) 
      Error("PowerMod: bad args");

   if (g.size() <= 1) {
      long gg = to_long(g);
      if (gg == 2)
         TwoPowerMod(h, e, F);
      else
         SPPowerMod(h, to_long(g), e, F);
      return;
   }

   if (e == 0) {
      set(h);
      return;
   }

   if (e == 1) {
      h = g;
      return;
   }

   if (e == -1) {
      InvMod(h, g, F);
      return;
   }

   if (e == 2) {
      SqrMod(h, g, F);
      return;
   }

   if (e == -2) {
      ZZ tmp;
      SqrMod(tmp, g, F);
      InvMod(h, tmp, F);
      return;
   }

   long n = NumBits(e);

   ZZ res;
   set(res);

   long i;

   if (n < 16) {
      // plain square-and-multiply algorithm

      for (i = n - 1; i >= 0; i--) {
         SqrMod(res, res, F);
         if (bit(e, i))
            MulMod(res, res, g, F);
      }

      if (e < 0) InvMod(res, res, F);

      h = res;
      return;
   }

   long k = OptWinSize(n);

   k = min(k, 5);

   vec_ZZ v;

   v.SetLength(1L << (k-1));

   v[0] = g;
 
   if (k > 1) {
      ZZ t;
      SqrMod(t, g, F);

      for (i = 1; i < (1L << (k-1)); i++)
         MulMod(v[i], v[i-1], t, F);
   }


   long val;
   long cnt;
   long m;

   val = 0;
   for (i = n-1; i >= 0; i--) {
      val = (val << 1) | bit(e, i); 
      if (val == 0)
         SqrMod(res, res, F);
      else if (val >= (1L << (k-1)) || i == 0) {
         cnt = 0;
         while ((val & 1) == 0) {
            val = val >> 1;
            cnt++;
         }

         m = val;
         while (m > 0) {
            SqrMod(res, res, F);
            m = m >> 1;
         }

         MulMod(res, res, v[val >> 1], F);

         while (cnt > 0) {
            SqrMod(res, res, F);
            cnt--;
         }

         val = 0;
      }
   }

   if (e < 0) InvMod(res, res, F);

   h = res;
}

void bit_and(ZZ& x, const ZZ& a, long b)
{
   ZZ B;
   conv(B, b);
   bit_and(x, a, B);
}

void bit_or(ZZ& x, const ZZ& a, long b)
{
   ZZ B;
   conv(B, b);
   bit_or(x, a, B);
}

void bit_xor(ZZ& x, const ZZ& a, long b)
{
   ZZ B;
   conv(B, b);
   bit_xor(x, a, B);
}


void power(ZZ& x, const ZZ& a, long e)
{
   if (a.size() > 1)
      zexp(a.rep, e, &x.rep);
   else
      zexps(to_long(a), e, &x.rep);
}


long power_long(long a, long e)
{
   if (e < 0) Error("power_long: negative exponent");

   if (e == 0) return 1;

   if (a == 1) return 1;
   if (a == -1) {
      if (e & 1)
         return -1;
      else
         return 1;
   }

   long res = 1;
   long i;

   for (i = 0; i < e; i++)
      res *= a;

   return res;
}

void AddMod(ZZ& x, const ZZ& a, const ZZ& b, const ZZ& n)
{
   if (&x != &n) 
      zaddmod(a.rep, b.rep, n.rep, &x.rep);
   else {
      ZZ tmp;
      zaddmod(a.rep, b.rep, n.rep, &tmp.rep);
      x = tmp;
   }
}


//  RANDOM NUMBER GENERATION

// Idea for this PRNG.  Iteratively hash seed using md5 
// to get 256 bytes to initialize arc4.
// Then take bytes from arc4 64 at a time, and feed these
// into md5, getting 4 words out.  These 4 words are placed into the output 
// stream and are also chained back into the md5.
// This should give high quality output, is not extremely fast,
// but fast enough.





// I make use of the md5 compression function,
// which I've modified to work on 64-bit machines


/*
 *  BEGIN RSA's md5 stuff
 *
 */

/*
 **********************************************************************
 ** md5.c                                                            **
 ** RSA Data Security, Inc. MD5 Message Digest Algorithm             **
 ** Created: 2/17/90 RLR                                             **
 ** Revised: 1/91 SRD,AJ,BSK,JT Reference C Version                  **
 **********************************************************************
 */

/*
 **********************************************************************
 ** Copyright (C) 1990, RSA Data Security, Inc. All rights reserved. **
 **                                                                  **
 ** License to copy and use this software is granted provided that   **
 ** it is identified as the "RSA Data Security, Inc. MD5 Message     **
 ** Digest Algorithm" in all material mentioning or referencing this **
 ** software or this function.                                       **
 **                                                                  **
 ** License is also granted to make and use derivative works         **
 ** provided that such works are identified as "derived from the RSA **
 ** Data Security, Inc. MD5 Message Digest Algorithm" in all         **
 ** material mentioning or referencing the derived work.             **
 **                                                                  **
 ** RSA Data Security, Inc. makes no representations concerning      **
 ** either the merchantability of this software or the suitability   **
 ** of this software for any particular purpose.  It is provided "as **
 ** is" without express or implied warranty of any kind.             **
 **                                                                  **
 ** These notices must be retained in any copies of any part of this **
 ** documentation and/or software.                                   **
 **********************************************************************
 */


#if (NTL_BITS_PER_LONG <= 32)
#define TRUNC32(x) (x)
#else
#define TRUNC32(x) ((x) & ((1UL << 32)-1UL))
#endif

/* F, G and H are basic MD5 functions: selection, majority, parity */
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) (TRUNC32((y) ^ ((x) | (~z)))) 

/* ROTATE_LEFT rotates x left n bits */
#define ROTATE_LEFT(x, n) (TRUNC32(((x) << (n)) | ((x) >> (32-(n)))))

/* FF, GG, HH, and II transformations for rounds 1, 2, 3, and 4 */
/* Rotation is separate from addition to prevent recomputation */
#define FF(a, b, c, d, x, s, ac) \
  {(a) = TRUNC32((a) + F((b), (c), (d)) + (x) + (ac)); \
   (a) = ROTATE_LEFT((a), (s)); \
   (a) = TRUNC32((a) + (b)); \
  }
#define GG(a, b, c, d, x, s, ac) \
  {(a) = TRUNC32((a) + G((b), (c), (d)) + (x) + (ac)); \
   (a) = ROTATE_LEFT((a), (s)); \
   (a) = TRUNC32((a) + (b)); \
  }
#define HH(a, b, c, d, x, s, ac) \
  {(a) = TRUNC32((a) + H((b), (c), (d)) + (x) + (ac)); \
   (a) = ROTATE_LEFT((a), (s)); \
   (a) = TRUNC32((a) + (b)); \
  }
#define II(a, b, c, d, x, s, ac) \
  {(a) = TRUNC32((a) + I((b), (c), (d)) + (x) + (ac)); \
   (a) = ROTATE_LEFT((a), (s)); \
   (a) = TRUNC32((a) + (b)); \
  }



static
void MD5_default_IV(unsigned long *buf)
{
   buf[0] = 0x67452301UL;
   buf[1] = 0xefcdab89UL;
   buf[2] = 0x98badcfeUL;
   buf[3] = 0x10325476UL;
}



/* Basic MD5 step. Transform buf based on in.
 */

static
void MD5_compress(unsigned long *buf, unsigned long *in)
{
  unsigned long a = buf[0], b = buf[1], c = buf[2], d = buf[3];

  /* Round 1 */
#define S11 7
#define S12 12
#define S13 17
#define S14 22
  FF ( a, b, c, d, in[ 0], S11, 3614090360UL); /* 1 */
  FF ( d, a, b, c, in[ 1], S12, 3905402710UL); /* 2 */
  FF ( c, d, a, b, in[ 2], S13,  606105819UL); /* 3 */
  FF ( b, c, d, a, in[ 3], S14, 3250441966UL); /* 4 */
  FF ( a, b, c, d, in[ 4], S11, 4118548399UL); /* 5 */
  FF ( d, a, b, c, in[ 5], S12, 1200080426UL); /* 6 */
  FF ( c, d, a, b, in[ 6], S13, 2821735955UL); /* 7 */
  FF ( b, c, d, a, in[ 7], S14, 4249261313UL); /* 8 */
  FF ( a, b, c, d, in[ 8], S11, 1770035416UL); /* 9 */
  FF ( d, a, b, c, in[ 9], S12, 2336552879UL); /* 10 */
  FF ( c, d, a, b, in[10], S13, 4294925233UL); /* 11 */
  FF ( b, c, d, a, in[11], S14, 2304563134UL); /* 12 */
  FF ( a, b, c, d, in[12], S11, 1804603682UL); /* 13 */
  FF ( d, a, b, c, in[13], S12, 4254626195UL); /* 14 */
  FF ( c, d, a, b, in[14], S13, 2792965006UL); /* 15 */
  FF ( b, c, d, a, in[15], S14, 1236535329UL); /* 16 */

  /* Round 2 */
#define S21 5
#define S22 9
#define S23 14
#define S24 20
  GG ( a, b, c, d, in[ 1], S21, 4129170786UL); /* 17 */
  GG ( d, a, b, c, in[ 6], S22, 3225465664UL); /* 18 */
  GG ( c, d, a, b, in[11], S23,  643717713UL); /* 19 */
  GG ( b, c, d, a, in[ 0], S24, 3921069994UL); /* 20 */
  GG ( a, b, c, d, in[ 5], S21, 3593408605UL); /* 21 */
  GG ( d, a, b, c, in[10], S22,   38016083UL); /* 22 */
  GG ( c, d, a, b, in[15], S23, 3634488961UL); /* 23 */
  GG ( b, c, d, a, in[ 4], S24, 3889429448UL); /* 24 */
  GG ( a, b, c, d, in[ 9], S21,  568446438UL); /* 25 */
  GG ( d, a, b, c, in[14], S22, 3275163606UL); /* 26 */
  GG ( c, d, a, b, in[ 3], S23, 4107603335UL); /* 27 */
  GG ( b, c, d, a, in[ 8], S24, 1163531501UL); /* 28 */
  GG ( a, b, c, d, in[13], S21, 2850285829UL); /* 29 */
  GG ( d, a, b, c, in[ 2], S22, 4243563512UL); /* 30 */
  GG ( c, d, a, b, in[ 7], S23, 1735328473UL); /* 31 */
  GG ( b, c, d, a, in[12], S24, 2368359562UL); /* 32 */

  /* Round 3 */
#define S31 4
#define S32 11
#define S33 16
#define S34 23
  HH ( a, b, c, d, in[ 5], S31, 4294588738UL); /* 33 */
  HH ( d, a, b, c, in[ 8], S32, 2272392833UL); /* 34 */
  HH ( c, d, a, b, in[11], S33, 1839030562UL); /* 35 */
  HH ( b, c, d, a, in[14], S34, 4259657740UL); /* 36 */
  HH ( a, b, c, d, in[ 1], S31, 2763975236UL); /* 37 */
  HH ( d, a, b, c, in[ 4], S32, 1272893353UL); /* 38 */
  HH ( c, d, a, b, in[ 7], S33, 4139469664UL); /* 39 */
  HH ( b, c, d, a, in[10], S34, 3200236656UL); /* 40 */
  HH ( a, b, c, d, in[13], S31,  681279174UL); /* 41 */
  HH ( d, a, b, c, in[ 0], S32, 3936430074UL); /* 42 */
  HH ( c, d, a, b, in[ 3], S33, 3572445317UL); /* 43 */
  HH ( b, c, d, a, in[ 6], S34,   76029189UL); /* 44 */
  HH ( a, b, c, d, in[ 9], S31, 3654602809UL); /* 45 */
  HH ( d, a, b, c, in[12], S32, 3873151461UL); /* 46 */
  HH ( c, d, a, b, in[15], S33,  530742520UL); /* 47 */
  HH ( b, c, d, a, in[ 2], S34, 3299628645UL); /* 48 */

  /* Round 4 */
#define S41 6
#define S42 10
#define S43 15
#define S44 21
  II ( a, b, c, d, in[ 0], S41, 4096336452UL); /* 49 */
  II ( d, a, b, c, in[ 7], S42, 1126891415UL); /* 50 */
  II ( c, d, a, b, in[14], S43, 2878612391UL); /* 51 */
  II ( b, c, d, a, in[ 5], S44, 4237533241UL); /* 52 */
  II ( a, b, c, d, in[12], S41, 1700485571UL); /* 53 */
  II ( d, a, b, c, in[ 3], S42, 2399980690UL); /* 54 */
  II ( c, d, a, b, in[10], S43, 4293915773UL); /* 55 */
  II ( b, c, d, a, in[ 1], S44, 2240044497UL); /* 56 */
  II ( a, b, c, d, in[ 8], S41, 1873313359UL); /* 57 */
  II ( d, a, b, c, in[15], S42, 4264355552UL); /* 58 */
  II ( c, d, a, b, in[ 6], S43, 2734768916UL); /* 59 */
  II ( b, c, d, a, in[13], S44, 1309151649UL); /* 60 */
  II ( a, b, c, d, in[ 4], S41, 4149444226UL); /* 61 */
  II ( d, a, b, c, in[11], S42, 3174756917UL); /* 62 */
  II ( c, d, a, b, in[ 2], S43,  718787259UL); /* 63 */
  II ( b, c, d, a, in[ 9], S44, 3951481745UL); /* 64 */

  buf[0] = TRUNC32(buf[0] + a);
  buf[1] = TRUNC32(buf[1] + b);
  buf[2] = TRUNC32(buf[2] + c);
  buf[3] = TRUNC32(buf[3] + d);
}


/*
 *  END RSA's md5 stuff
 *
 */


static
void words_from_bytes(unsigned long *txtl, unsigned char *txtc, long n)
{
   long i;
   unsigned long v;

   for (i = 0; i < n; i++) {
      v = txtc[4*i];
      v += ((unsigned long) (txtc[4*i+1])) << 8;
      v += ((unsigned long) (txtc[4*i+2])) << 16;
      v += ((unsigned long) (txtc[4*i+3])) << 24;
      txtl[i] = v;
   }
}

static 
void bytes_from_words(unsigned char *txtc, unsigned long *txtl, long n)
{
   long i;
   unsigned long v;

   for (i = 0; i < n; i++) {
      v = txtl[i];
      txtc[4*i] = v & 255;
      v = v >> 8;
      txtc[4*i+1] = v & 255;
      v = v >> 8;
      txtc[4*i+2] = v & 255;
      v = v >> 8;
      txtc[4*i+3] = v & 255;
   }
}


static
void MD5_compress1(unsigned long *buf, unsigned char *in, long n)
{
   unsigned long txtl[16];
   unsigned char txtc[64]; 
   long i, j, k;

   if (n < 0) n = 0;

   i = 0;
   while (i < n) {
      k = n-i;
      if (k > 64) k = 64;
      for (j = 0; j < k; j++)
         txtc[j] = in[i+j];
      for (; j < 64; j++)
         txtc[j] = 0;
      words_from_bytes(txtl, txtc, 16);
      MD5_compress(buf, txtl);
      i += k;
   }
}


// the "cipherpunk" version of arc4 

struct arc4_key
{      
    unsigned char state[256];       
    unsigned char x;        
    unsigned char y;
};


inline
void swap_byte(unsigned char *a, unsigned char *b)
{
    unsigned char swapByte; 
    
    swapByte = *a; 
    *a = *b;      
    *b = swapByte;
}

static
void prepare_key(unsigned char *key_data_ptr, 
                 long key_data_len, arc4_key *key)
{
    unsigned char index1;
    unsigned char index2;
    unsigned char* state;
    long counter;     
    
    state = &key->state[0];         
    for(counter = 0; counter < 256; counter++)              
       state[counter] = counter;               
    key->x = 0;     
    key->y = 0;     
    index1 = 0;     
    index2 = 0;             
    for(counter = 0; counter < 256; counter++)      
    {               
         index2 = (key_data_ptr[index1] + state[counter] + index2) & 255;                
         swap_byte(&state[counter], &state[index2]);            

         index1 = (index1 + 1) % key_data_len;  
    }       
}



static
void arc4(unsigned char *buffer_ptr, long buffer_len, arc4_key *key)
{ 
    unsigned char x;
    unsigned char y;
    unsigned char* state;
    unsigned char xorIndex;
    long counter;              
    
    x = key->x;     
    y = key->y;     
    
    state = &key->state[0];         
    for(counter = 0; counter < buffer_len; counter ++)      
    {               
         x = (x + 1) & 255;
         y = (state[x] + y) & 255;
         swap_byte(&state[x], &state[y]);                        
              
         xorIndex = (state[x] + state[y]) & 255;
              
         buffer_ptr[counter] = state[xorIndex];         
     }               
     key->x = x;     
     key->y = y;
}

// global state information for PRNG

static long ran_initialized = 0;
static unsigned long ran_buf[4];
static long ran_count;
static arc4_key ran_key;

static unsigned long default_md5_tab[16] = {
744663023UL, 1011602954UL, 3163087192UL, 3383838527UL, 
3305324122UL, 3197458079UL, 2266495600UL, 2760303563UL, 
346234297UL, 1919920720UL, 1896169861UL, 2192176675UL, 
2027150322UL, 2090160759UL, 2134858730UL, 1131796244UL
};


void build_arc4_tab(unsigned char *seed_bytes, const ZZ& s)
{
   long nb = NumBytes(s);
   
   unsigned char *txt;

   typedef unsigned char u_char;
   txt = new u_char[nb + 64];
   if (!txt) Error("out of memory");

   BytesFromZZ(txt, s, nb);

   bytes_from_words(txt+nb, default_md5_tab, 16);

   unsigned long buf[4];
   MD5_default_IV(buf);

   long i;
   for (i = 0; i < 16; i++) {
      MD5_compress1(buf, txt, nb + 64);
      bytes_from_words(seed_bytes + 16*i, buf, 4);
   }

   delete [] txt;
}


void SetSeed(const ZZ& s)
{
   unsigned char seed_bytes[256];
   build_arc4_tab(seed_bytes, s);
   prepare_key(seed_bytes, 256, &ran_key);

   MD5_default_IV(ran_buf);

   ran_initialized = 1;
   ran_count = 4;
}

static 
unsigned long ran_next_32bits()
{
   unsigned long wbuf[16];
   unsigned char bbuf[64];

   unsigned long res;

   if (!ran_initialized) SetSeed(ZZ::zero());

   if (ran_count >= 4) {
      arc4(bbuf, 64, &ran_key);
      words_from_bytes(wbuf, bbuf, 16);
      MD5_compress(ran_buf, wbuf);
      ran_count = 0;
   }

   res = ran_buf[ran_count];
   ran_count++;
   return res;
}




#if (NTL_BITS_PER_LONG <= 32)


unsigned long RandomWord()
{
   return ran_next_32bits();
}

#else


unsigned long RandomWord()
{
   const long n = (NTL_BITS_PER_LONG+31)/32;
   long i;
   unsigned long res = 0;
   for (i = 0; i < n; i++) 
      res = (res << 32) ^ ran_next_32bits();

   return res;
}

#endif





long RandomLen_long(long l)
{
   if (l <= 0) return 0;
   if (l == 1) return 1;
   if (l >= NTL_BITS_PER_LONG) 
      Error("RandomLen: l out out of range");
   return long(RandomWord() & ((1UL << (l-1))-1UL)) + (1L << (l-1)); 
}


long RandomBits_long(long l)
{
   if (l <= 0) return 0;
   if (l >= NTL_BITS_PER_LONG) 
      Error("RandomLen: l out out of range");
   return long(RandomWord() & ((1UL << l)-1UL));
}


void RandomBits(ZZ& x, long bitlength)
{
   verylong *aa = &x.rep;
   long wc, bc, i;
   verylong a;

   if (bitlength <= 0) {
      zzero(aa);
      return;
   }

   wc = bitlength /NTL_NBITS;
   bc = bitlength - wc*NTL_NBITS;

   if (bc != 0) 
      wc++;
   else
      bc = NTL_NBITS;

   a = *aa;

   zsetlength(&a, wc);
   *aa = a;

   for (i = 1; i <= wc; i++)
      a[i] = long(RandomWord() & NTL_RADIXM);

   a[wc] &= ((1L << bc)-1L);

   while (wc > 1 && a[wc] == 0) wc--;
   a[0] = wc;
}

void RandomLen(ZZ& x, long bitlength)
{
   verylong *aa = &x.rep;
   long wc, bc, i;
   verylong a;

   if (bitlength <= 0) {
      zzero(aa);
      return;
   }

   if (bitlength == 1) {
      zone(aa);
      return;
   }

   wc = bitlength /NTL_NBITS;
   bc = bitlength - wc*NTL_NBITS;

   if (bc != 0) 
      wc++;
   else
      bc = NTL_NBITS;

   a = *aa;

   zsetlength(&a, wc);
   *aa = a;

   a[0] = wc;
   for (i = 1; i <= wc; i++)
      a[i] = long(RandomWord() & NTL_RADIXM);

   a[wc] = (((1L << (bc-1))-1) & a[wc]) | (1L << (bc-1));
}

void RandomBnd(ZZ& x, const ZZ& bnd)
{
   if (bnd <= 1) {
      x = 0;
      return;
   }

   long k = NumBits(bnd);

   if (weight(bnd) == 1) {
      RandomBits(x, k-1);
      return;
   }

   long l = k + 8;

   ZZ t, r, t1;

   do {
      RandomBits(t, l);
      r = t;
      QuickRem(r, bnd);
      SubPos(t1, bnd, r);
      add(t, t, t1);
   } while (NumBits(t) > l);

   x = r;
}

long RandomBnd(long bnd)
{
   if (bnd <= 1) return 0;

   long k = NumBits(bnd);
   long l = k + 4;

   if (l > NTL_BITS_PER_LONG-2) {
      ZZ Bnd, res;

      Bnd = bnd;
      RandomBnd(res, Bnd);
      return to_long(res);
   }

   long t, r;

   do {
      t = RandomBits_long(l);
      r = t % bnd;
   } while (t + bnd - r > (1L << l)); 

   return r;
}



void ZZFromBytes(ZZ& x, const unsigned char *p, long n)
{
   if (n <= 0) {
      x = 0;
      return;
   }

   if (n > (NTL_MAX_LONG-(NTL_NBITS-1))/8)
      Error("ZZFromBytes: excessive length");

   long sz = (n*8 + NTL_NBITS-1)/NTL_NBITS;

   x.SetSize(sz);

   long i;

   verylong a = x.rep;

   for (i = 1; i <= sz; i++)
      a[i] = 0;

   for (i = 0; i < n; i++) {
      long bitpos = i*8;
      long wordpos = bitpos/NTL_NBITS;
      long bitoffset = bitpos - wordpos*NTL_NBITS;

      a[wordpos+1] |= ((long(p[i]) & 255) << bitoffset);

      long diff = NTL_NBITS-bitoffset;

      if (diff < 8) {
         a[wordpos+1] &= NTL_RADIXM;
         a[wordpos+2] = (long(p[i]) & 255) >> diff;
      }
   }

   while (sz > 1 && a[sz] == 0) sz--;
   a[0] = sz;
}

void BytesFromZZ(unsigned char *p, const ZZ& s, long nn)
{
   long k = NumBits(s);
   long n = (k+7)/8;

   long sz = s.size();

   verylong a = s.rep;
 
   long i;

   long min_n = min(n, nn);

   for (i = 0; i < min_n; i++) {
      long bitpos = i*8;
      long wordpos = bitpos/NTL_NBITS;
      long bitoffset = bitpos - wordpos*NTL_NBITS;

      p[i] = (a[wordpos+1] >> bitoffset) & 255;

      long diff = NTL_NBITS - bitoffset;

      if (diff < 8 && wordpos < sz-1) {
         long msk = (1L << (8-diff))-1;
         p[i] |= ((a[wordpos+2] & msk) << diff);
      }
   }

   for (i = min_n; i < nn; i++)
      p[i] = 0;
}
