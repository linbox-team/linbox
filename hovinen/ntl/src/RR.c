
#include <NTL/RR.h>
#include <NTL/IsFinite.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>


long RR::prec = 150;

void RR::SetPrecision(long p)
{
   if (p < 53)
      p = 53;

   if (p >= (1L << (NTL_BITS_PER_LONG-4)))
      Error("RR: precision too high");

   prec = p;
}

long RR::oprec = 10;

void RR::SetOutputPrecision(long p)
{
   if (p < 1)
      p = 1;

   oprec = p;
}



void normalize(RR& z, const RR& y, long residual = 0)
{
   long len = NumBits(y.x);

   if (len > RR::prec) {
      const long *a = y.x.rep;

      long sgn;

      long direction;
      long p, wh, bl;

      if (a[0] > 0)
         sgn = 1;
      else
         sgn = -1;

      p = len - RR::prec - 1;
      bl = (p/NTL_NBITS);
      wh = 1L << (p - NTL_NBITS*bl);
      bl++;

      if (a[bl] & wh) {
         // bit is 1...we have to see if lower bits are all 0
         // in order to implement "round to even"

         if (a[bl] & (wh - 1)) 
            direction = 1;
         else {
            long i = bl - 1;
            while (i > 0 && a[i] == 0) i--;
            if (i > 0)
               direction = 1;
            else
               direction = 0;
         }

         // use residual to break ties

         if (direction == 0 && residual != 0) {
            if (residual == sgn)
               direction = 1;
            else 
               direction = -1;
         }

         if (direction == 0) {
            // round to even

            wh = wh << 1;
            if (wh == NTL_RADIX) {
               wh = 1;
               bl++;
            }

            if (a[bl] & wh)
               direction = 1;
            else
               direction = -1;
         }
      }
      else
         direction = -1;

      RightShift(z.x, y.x, len - RR::prec);
      if (direction == 1) 
         add(z.x, z.x, sgn);

      z.e = y.e + len - RR::prec;
   }
   else if (len == 0) {
      clear(z.x);
      z.e = 0;
   }
   else {
      z.x = y.x;
      z.e = y.e;
   }

   if (!IsOdd(z.x))
      z.e += MakeOdd(z.x);

   if (z.e >= (1L << (NTL_BITS_PER_LONG-4)))
      Error("RR: overflow");

   if (z.e <= -(1L << (NTL_BITS_PER_LONG-4)))
      Error("RR: underflow");
}

inline void xcopy(RR& x, const RR& a)
   { normalize(x, a); }

// xcopy emulates old assignment semantics...
// many routines here implicitly assume assignment normalizes,
// but that is no longer the case as of v3.0.

void conv(RR& x, const RR& a)
{
   normalize(x, a);
}


long IsZero(const RR& a)
{
   return IsZero(a.x);
}

long IsOne(const RR& a)
{
   return a.e == 0 && IsOne(a.x);
}

long sign(const RR& a)
{
   return sign(a.x);
}

void clear(RR& z)
{
   z.e = 0;
   clear(z.x);
}

void set(RR& z)
{
   z.e = 0;
   set(z.x);
}


void add(RR& z, const RR& a, const RR& b)
{
   static RR t;

   if (IsZero(a.x)) {
      xcopy(z, b);
      return;
   }

   if (IsZero(b.x)) {
      xcopy(z, a);
      return;
   }

   if (a.e > b.e) {
      if (a.e-b.e - max(RR::prec-NumBits(a.x),0) >= NumBits(b.x) + 2)
         normalize(z, a, sign(b));
      else {
         LeftShift(t.x, a.x, a.e-b.e);
         add(t.x, t.x, b.x);
         t.e = b.e;
         normalize(z, t);
      }
   }
   else if (a.e < b.e) {
      if (b.e-a.e - max(RR::prec-NumBits(b.x),0) >= NumBits(a.x) + 2)
         normalize(z, b, sign(a));
      else {
         LeftShift(t.x, b.x, b.e-a.e);
         add(t.x, t.x, a.x);
         t.e = a.e;
         normalize(z, t);
      }
   }
   else {
      add(t.x, a.x, b.x);
      t.e = a.e;
      normalize(z, t);
   }
}

void sub(RR& z, const RR& a, const RR& b)
{
   static RR t;

   if (IsZero(a.x)) {
      negate(z, b);
      return;
   }

   if (IsZero(b.x)) {
      xcopy(z, a);
      return;
   }

   if (a.e > b.e) {
      if (a.e-b.e - max(RR::prec-NumBits(a.x),0) >= NumBits(b.x) + 2)
         normalize(z, a, -sign(b));
      else {
         LeftShift(t.x, a.x, a.e-b.e);
         sub(t.x, t.x, b.x);
         t.e = b.e;
         xcopy(z, t);
      }
   }
   else if (a.e < b.e) {
      if (b.e-a.e - max(RR::prec-NumBits(b.x),0) >= NumBits(a.x) + 2) {
         normalize(z, b, -sign(a));
         negate(z.x, z.x);
      }
      else {
         LeftShift(t.x, b.x, b.e-a.e);
         sub(t.x, a.x, t.x);
         t.e = a.e;
         xcopy(z, t);
      }
   }
   else {
      sub(t.x, a.x, b.x);
      t.e = a.e;
      normalize(z, t);
   }
}

void negate(RR& z, const RR& a)
{
   xcopy(z, a);
   negate(z.x, z.x);
}

void abs(RR& z, const RR& a)
{
   xcopy(z, a);
   abs(z.x, z.x);
}


void mul(RR& z, const RR& a, const RR& b)
{
   static RR t;

   mul(t.x, a.x, b.x);
   t.e = a.e + b.e;
   xcopy(z, t);
}

void sqr(RR& z, const RR& a)
{
   static RR t;

   sqr(t.x, a.x);
   t.e = a.e + a.e;
   xcopy(z, t);
}

void div(RR& z, const RR& a, const RR& b)
{
   if (IsZero(b))
      Error("RR: division by zero");

   if (IsZero(a)) {
      clear(z);
      return;
   }

   long la = NumBits(a.x);
   long lb = NumBits(b.x);

   long neg = (sign(a) != sign(b));

   long k = RR::prec - la + lb + 1;
   if (k < 0) k = 0;

   static RR t;
   static ZZ A, B, R;

   abs(A, a.x);
   LeftShift(A, A, k);

   abs(B, b.x);
   DivRem(t.x, R, A, B);

   t.e = a.e - b.e - k;

   normalize(z, t, !IsZero(R));

   if (neg)
      negate(z.x, z.x);
}

void SqrRoot(RR& z, const RR& a)
{
   if (sign(a) < 0)
      Error("RR: attempt to take square root of negative number");

   if (IsZero(a)) {
      clear(z);
      return;
   }

   RR t;
   ZZ T1, T2;
   long k;

   k = 2*RR::prec - NumBits(a.x) + 1;

   if (k < 0) k = 0;

   if ((a.e - k) & 1) k++;

   LeftShift(T1, a.x, k);
   SqrRoot(t.x, T1);
   t.e = (a.e - k)/2;
   sqr(T2, T1);

   normalize(z, t, T2 < T1);
}
   




void swap(RR& a, RR& b)
{
   swap(a.x, b.x);
   swap(a.e, b.e);
}

long compare(const RR& a, const RR& b)
{
   static RR t;

   sub(t, a, b);
   return sign(t);
}



long operator==(const RR& a, const RR& b) 
{
   return a.e == b.e && a.x == b.x;
}


void trunc(RR& z, const RR& a)
{
   static RR t;

   if (a.e >= 0) 
      xcopy(z, a);
   else {
      RightShift(t.x, a.x, -a.e);
      t.e = 0;
      xcopy(z, t);
   }
}

void floor(RR& z, const RR& a)
{
   static RR t;

   if (a.e >= 0) 
      xcopy(z, a);
   else {
      RightShift(t.x, a.x, -a.e);
      if (sign(a.x) < 0)
         add(t.x, t.x, -1);
      t.e = 0;
      xcopy(z, t);
   }
}

void ceil(RR& z, const RR& a)
{
   static RR t;

   if (a.e >= 0)
      xcopy(z, a);
   else {
      RightShift(t.x, a.x, -a.e);
      if (sign(a.x) > 0)
         add(t.x, t.x, 1);
      t.e = 0;
      xcopy(z, t);
   }
}
   

void conv(RR& z, const ZZ& a)
{
   static RR t;

   t.x = a;
   t.e = 0;

   xcopy(z, t);
}

void conv(RR& z, long a)
{
   if (a == 0) {
      clear(z);
      return;
   }

   if (a == 1) {
      set(z);
      return;
   }

   static ZZ t;
   t = a;
   conv(z, t);
}


void conv(RR& z, double a)
{
   if (a == 0) {
      clear(z);
      return;
   }

   if (a == 1) {
      set(z);
      return;
   }

   if (!IsFinite(&a))
      Error("RR: conversion of a non-finite double");

   int e;
   double f;
   static RR t;

   f = frexp(a, &e);

   f = f * NTL_FDOUBLE_PRECISION;
   f = f * 4;

   conv(t.x, f);
   t.e = e - (NTL_DOUBLE_PRECISION + 1);

   xcopy(z, t);
}


void conv(ZZ& z, const RR& a)
{
   static RR t;
   long old_p = RR::precision();
   RR::SetPrecision(NumBits(a.x) + 1); // prevents rounding

   floor(t, a);
   LeftShift(z, t.x, t.e);

   RR::SetPrecision(old_p);
}

void conv(long& z, const RR& a)
{
   static ZZ t;
   conv(t, a);
   conv(z, t);
}

void conv(double& z, const RR& aa)
{
   long old_p;
   double x;
   int e;
   static RR a;

   old_p = RR::prec;
   RR::prec = NTL_DOUBLE_PRECISION;

   xcopy(a, aa);  // round to NTL_DOUBLE_PRECISION bits to avoid double overflow

   e = a.e;
   if (e != a.e) Error("RR: overflow in conversion to double");
   conv(x, a.x);


   z = ldexp(x, e);

   RR::prec = old_p;
}


void add(RR& z, const RR& a, double b)
{
   static RR B;
   B = b;
   add(z, a, B);
}



void sub(RR& z, const RR& a, double b)
{
   static RR B;
   B = b;
   sub(z, a, B);
}

void sub(RR& z, double a, const RR& b)
{
   static RR A;
   A = a;
   sub(z, A, b);
}



void mul(RR& z, const RR& a, double b)
{
   static RR B;
   B = b;
   mul(z, a, B);
}


void div(RR& z, const RR& a, double b)
{
   static RR B;
   B = b;
   div(z, a, B);
}

void div(RR& z, double a, const RR& b)
{
   static RR A;
   A = a;
   div(z, A, b);
}


void inv(RR& z, const RR& a)
{
   static RR one = to_RR(1);
   div(z, one, a);
}


long compare(const RR& a, double b)
{
   if (b == 0) return sign(a);

   static RR B;
   B = b;
   return compare(a, B);
}


long operator==(const RR& a, double b) 
{
   if (b == 0) return IsZero(a);
   if (b == 1) return IsOne(a);

   static RR B;
   B = b;
   return a == B;
}


void power(RR& z, const RR& a, long e)
{
   RR b, res;
   long neg;

   long n = NumBits(e);

   long p = RR::precision();
   RR::SetPrecision(p + n + 10);

   xcopy(b, a);
   if (e < 0) {
      e = -e;
      neg = 1;
   }
   else
      neg = 0;

   set(res);
   long i;

   for (i = n-1; i >= 0; i--) {
      sqr(res, res);
      if (bit(e, i))
         mul(res, res, b);
   }

   if (neg) 
      inv(z, res);
   else
      xcopy(z, res);

   RR::SetPrecision(p);
}

ostream& operator<<(ostream& s, const RR& a)
{
   if (IsZero(a)) {
      s << "0";
      return s;
   }

   long old_p = RR::precision();

   RR::SetPrecision(long(RR::OutputPrecision()*3.321928095) + 20);

   RR b;
   long neg;

   if (a < 0) {
      negate(b, a);
      neg = 1;
   }
   else {
      xcopy(b, a);
      neg = 0;
   }

   long k;

   k = long((RR::OutputPrecision()*3.321928095-NumBits(b.mantissa())
            -b.exponent()) / 3.321928095);


   RR c;

   power(c, to_RR(10), k);
   mul(b, b, c);

   power(c, to_RR(10), RR::OutputPrecision());

   while (b < c) {
      mul(b, b, 10);
      k++;
   }

   while (b >= c) {
      div(b, b, 10);
      k--;
   }

   add(b, b, 0.5);
   k = -k;

   ZZ B;
   conv(B, b);

   char *bp = new char[RR::OutputPrecision()+10];

   if (!bp) Error("RR output: out of memory");

   long len, i;

   len = 0;
   do {
      bp[len] = DivRem(B, B, 10) + '0';
      len++;
   } while (B > 0);

   for (i = 0; i < len/2; i++) {
      char tmp;
      tmp = bp[i];
      bp[i] = bp[len-1-i];
      bp[len-1-i] = tmp;
   }

   i = len-1;
   while (bp[i] == '0') i--;

   k += (len-1-i);
   len = i+1;

   bp[len] = '\0';

   if (k > 3 || k < -len - 3) {
      // use scientific notation

      if (neg) s << "-";
      s << "0." << bp << "e" << (k + len);
   }
   else if (k >= 0) {
      if (neg) s << "-";
      s << bp;
      for (i = 0; i < k; i++) 
         s << "0";
   }
   else if (k <= -len) {
      if (neg) s << "-";
      s << "0.";
      for (i = 0; i < -len-k; i++)
         s << "0";
      s << bp;
   }
   else {
      if (neg) s << "-";
      for (i = 0; i < len+k; i++)
         s << bp[i];

      s << ".";

      for (i = len+k; i < len; i++)
         s << bp[i];
   }

   RR::SetPrecision(old_p);
   delete [] bp;
   return s;
}

istream& operator>>(istream& s, RR& x)
{
   long c;
   long sign;
   ZZ a, b;

   if (!s) Error("bad RR input");


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

   long got1 = 0;
   long got_dot = 0;
   long got2 = 0;

   a = 0;
   b = 1;

   if (c >= '0' && c <= '9') {
      got1 = 1;

      while (c >= '0' && c <= '9') {
         mul(a, a, 10);
         add(a, a, c-'0');
         s.get();
         c = s.peek();
      }
   }

   if (c == '.') {
      got_dot = 1;

      s.get();
      c = s.peek();

      if (c >= '0' && c <= '9') {
         got2 = 1;
   
         while (c >= '0' && c <= '9') {
            mul(a, a, 10);
            add(a, a, c-'0');
            mul(b, b, 10);
            s.get();
            c = s.peek();
         }
      }
   }

   if (got_dot && !got1 && !got2)  Error("bad RR input");

   ZZ e;

   long got_e = 0;
   long e_sign;

   if (c == 'e' || c == 'E') {
      got_e = 1;

      s.get();
      c = s.peek();

      if (c == '-') {
         e_sign = -1;
         s.get();
         c = s.peek();
      }
      else if (c == '+') {
         e_sign = 1;
         s.get();
         c = s.peek();
      }
      else
         e_sign = 1;

      if (c < '0' || c > '9') Error("bad RR input");

      e = 0;
      while (c >= '0' && c <= '9') {
         mul(e, e, 10);
         add(e, e, c-'0');
         s.get();
         c = s.peek();
      }
   }

   if (!got1 && !got2 && !got_e) Error("bad RR input");

   RR t1, t2, v;

   long old_p = RR::precision();

   if (got1 || got2) {
      RR::SetPrecision(max(NumBits(a), NumBits(b)));
      conv(t1, a);
      conv(t2, b);
      if (got_e)
         RR::SetPrecision(old_p + 10);
      else
         RR::SetPrecision(old_p);
      div(v, t1, t2);
   }
   else
      set(v);

   if (sign < 0)
      negate(v, v);

   if (got_e) {
      if (e >= (1L << (NTL_BITS_PER_LONG-4))) Error("RR input overflow");
      long E;
      conv(E, e);
      if (e_sign < 0) E = -E;
      RR::SetPrecision(old_p + 10);
      power(t1, to_RR(10), E);
      mul(v, v, t1);
      RR::SetPrecision(old_p);
   }

   xcopy(x, v);
   return s;
}


void conv(RR& z, const xdouble& a)
{
   static RR t1, t2;


   conv(t1, a.mantissa());

   if (a.exponent() >= 0)
      power(t2, to_RR(NTL_XD_BOUND), a.exponent());
   else
      power(t2, to_RR(1/NTL_XD_BOUND), -a.exponent());

   mul(z, t1, t2);
}


void conv(xdouble& z, const RR& a)
{
   xdouble x;
   xdouble y;

   conv(x, a.x);
   power2(y, a.e);
   z = x*y;
}
      
void power2(RR& z, long e)
{
   if (e >= (1L << (NTL_BITS_PER_LONG-4)))
      Error("RR: overflow");

   if (e <= -(1L << (NTL_BITS_PER_LONG-4)))
      Error("RR: underflow");

   set(z.x); 
   z.e = e;
}

void conv(RR& z, const quad_float& a)
{
   static RR hi, lo, res;

   long old_p = RR::prec;
   RR::prec = 2*NTL_DOUBLE_PRECISION;

   conv(hi, a.hi);
   conv(lo, a.lo);

   add(res, hi, lo);

   RR::prec = old_p;

   xcopy(z, res);
}


void conv(quad_float& z, const RR& aa)
{
   long old_p;
   quad_float x, y;
   static RR a;

   old_p = RR::prec;
   RR::prec = 2*NTL_DOUBLE_PRECISION;

   xcopy(a, aa);  // round to 2*NTL_DOUBLE_PRECISION bits to avoid double overflow

   conv(x, a.x);
   power2(y, a.e);

   z = x*y;

   RR::prec = old_p;
}

void conv(RR& x, const char *s)
{
   long c;
   long sign;
   ZZ a, b;
   long i = 0;

   if (!s) Error("bad RR input");


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

   long got1 = 0;
   long got_dot = 0;
   long got2 = 0;

   a = 0;
   b = 1;

   if (c >= '0' && c <= '9') {
      got1 = 1;

      while (c >= '0' && c <= '9') {
         mul(a, a, 10);
         add(a, a, c-'0');
         i++;
         c = s[i];
      }
   }

   if (c == '.') {
      got_dot = 1;

      i++;
      c = s[i];

      if (c >= '0' && c <= '9') {
         got2 = 1;
   
         while (c >= '0' && c <= '9') {
            mul(a, a, 10);
            add(a, a, c-'0');
            mul(b, b, 10);
            i++;
            c = s[i];
         }
      }
   }

   if (got_dot && !got1 && !got2)  Error("bad RR input");

   ZZ e;

   long got_e = 0;
   long e_sign;

   if (c == 'e' || c == 'E') {
      got_e = 1;

      i++;
      c = s[i];

      if (c == '-') {
         e_sign = -1;
         i++;
         c = s[i];
      }
      else if (c == '+') {
         e_sign = 1;
         i++;
         c = s[i];
      }
      else
         e_sign = 1;

      if (c < '0' || c > '9') Error("bad RR input");

      e = 0;
      while (c >= '0' && c <= '9') {
         mul(e, e, 10);
         add(e, e, c-'0');
         i++;
         c = s[i];
      }
   }

   if (!got1 && !got2 && !got_e) Error("bad RR input");

   RR t1, t2, v;

   long old_p = RR::precision();

   if (got1 || got2) {
      RR::SetPrecision(max(NumBits(a), NumBits(b)));
      conv(t1, a);
      conv(t2, b);
      if (got_e)
         RR::SetPrecision(old_p + 10);
      else
         RR::SetPrecision(old_p);
      div(v, t1, t2);
   }
   else
      set(v);

   if (sign < 0)
      negate(v, v);

   if (got_e) {
      if (e >= (1L << (NTL_BITS_PER_LONG-4))) Error("RR input overflow");
      long E;
      conv(E, e);
      if (e_sign < 0) E = -E;
      RR::SetPrecision(old_p + 10);
      power(t1, to_RR(10), E);
      mul(v, v, t1);
      RR::SetPrecision(old_p);
   }

   xcopy(x, v);
}


void ReallyComputeE(RR& res)
{
   long p = RR::precision();
   RR::SetPrecision(p + NumBits(p) + 10);

   RR s, s1, t;

   s = 1;
   t = 1;

   long i;

   for (i = 2; ; i++) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      div(t, t, i);
   }

   RR::SetPrecision(p);
   xcopy(res, s);
}

void ComputeE(RR& res)
{
   static long prec = 0;
   static RR e;

   long p = RR::precision();

   if (prec <= p + 10) {
      prec = p + 20;
      RR::SetPrecision(prec);
      ReallyComputeE(e);
      RR::SetPrecision(p);
   }

   xcopy(res, e);
}


void exp(RR& res, const RR& x)
{
   if (x >= (1L << (NTL_BITS_PER_LONG-4)) || x <= -(1L << (NTL_BITS_PER_LONG-4)))
      Error("RR: overflow");

   long p = RR::precision();

   // step 0: write x = n + f, n an integer and |f| <= 1/2
   //    careful -- we want to compute f to > p bits of precision


   RR::SetPrecision(p + 10);

   long n = to_long(x);
   RR f;
   sub(f, x, n);
   if (f > 0.5) {
      n++;
      sub(f, x, n);
   }

   // step 1: calculate t1 = e^n by repeated squaring

   RR::SetPrecision(p + NumBits(n) + 10);

   RR e, t1;
   ComputeE(e);
   power(t1, e, n); 

   // step 2: calculate t2 = e^f using Taylor series expansion

   RR::SetPrecision(p + NumBits(p) + 10);

   RR t2, s, s1, t;
   long i;

   s = 0;
   t = 1;

   for (i = 1; ; i++) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t, t, f);
      div(t, t, i);
   }

   xcopy(t2, s);

   RR::SetPrecision(p);

   mul(res, t1, t2);
}



void ReallyComputeLn2(RR& res)
{
   long p = RR::precision();
   RR::SetPrecision(p + NumBits(p) + 10);

   RR s, s1, t, t1;

   s = 0;
   t = 0.5;
   t1 = 0.5;

   long i;

   for (i = 2; ; i++) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t1, t1, 0.5);
      div(t, t1, i);
   }

   RR::SetPrecision(p);
   xcopy(res, s);
}


void ComputeLn2(RR& res)
{
   static long prec = 0;
   static RR ln2;

   long p = RR::precision();

   if (prec <= p + 10) {
      prec = p + 20;
      RR::SetPrecision(prec);
      ReallyComputeLn2(ln2);
      RR::SetPrecision(p);
   }

   xcopy(res, ln2);
}

long Lg2(const RR& x)
{
   return NumBits(x.mantissa()) + x.exponent();
}

void log(RR& res, const RR& x)
{
   if (x <= 0) Error("argument to log must be positive");

   long p = RR::precision();

   RR::SetPrecision(p + NumBits(p) + 10);

   RR y;
   long n;

   // re-write x = 2^n * (1 - y), where -1/2 < y < 1/4  (so 3/4 < 1-y < 3/2)

   if (x > 0.75 && x < 1.5) {
      n = 0;
      sub(y, 1, x);
   }
   else {
      n = Lg2(x) - 1;
      RR t;
      power2(t, -n);
      mul(t, t, x);
      while (t > 1.5) {
         mul(t, t, 0.5);
         n++;
      }

      sub(y, 1, t);
   }

   // compute s = - ln(1-y) by power series expansion

   RR s, s1, t, t1;

   s = 0;
   xcopy(t, y);
   xcopy(t1, y);

   long i;

   for (i = 2; ; i++) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t1, t1, y);
      div(t, t1, i);
   }

   if (n == 0) 
      t = 0;
   else {
      ComputeLn2(t);
      mul(t, t, n);
   }

   RR::SetPrecision(p);

   sub(res, t, s);
}


void ComputeLn10(RR& res)
{
   static long prec = 0;
   static RR ln10;

   long p = RR::precision();

   if (prec <= p + 10) {
      prec = p + 20;
      RR::SetPrecision(prec);
      log(ln10, to_RR(10));
      RR::SetPrecision(p);
   }

   xcopy(res, ln10);
}

void log10(RR& res, const RR& x)
{
   long p = RR::precision();
   RR::SetPrecision(p + 10);

   RR ln10, t1, t2;
   ComputeLn10(ln10);

   log(t1, x);
   div(t2, t1, ln10);

   RR::SetPrecision(p);

   xcopy(res, t2);
}


void expm1(RR& res, const RR& x)
{
   if (x < -0.5 || x > 0.5) {
      RR t;
      exp(t, x);
      sub(res, t, 1);
      return;
   }

   long p = RR::precision();

   RR::SetPrecision(p + NumBits(p) + 10);

   RR f;

   xcopy(f, x);

   RR s, s1, t;
   long i;

   s = 0;
   xcopy(t, f);

   for (i = 2; ; i++) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t, t, f);
      div(t, t, i);
   }

   RR::SetPrecision(p);

   xcopy(res, s);
}



void log1p(RR& res, const RR& x)
{
   if (x < -0.5 || x > 0.5) {
      log(res, 1 + x);
      return;
   }

   long p = RR::precision();

   RR::SetPrecision(p + NumBits(p) + 10);

   RR y;

   negate(y, x);

   // compute s = - ln(1-y) by power series expansion

   RR s, s1, t, t1;

   s = 0;
   xcopy(t, y);
   xcopy(t1, y);

   long i;

   for (i = 2; ; i++) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t1, t1, y);
      div(t, t1, i);
   }

   RR::SetPrecision(p);

   negate(res, s);

}


void pow(RR& res, const RR& x, const RR& y)
{
   if (y == 0) {
      res = 1;
      return;
   }

   if (x == 0) {
      res = 0;
      return;
   }

   if (x == 1) {
      res = 1;
      return;
   }

   if (x < 0) {
      Error("pow: sorry...first argument to pow must be nonnegative");
   }

   long p = RR::precision();

   // calculate working precison...one could use p + NTL_BITS_PER_LONG + 10,
   // for example, but we want the behaviour to be machine independent.
   // so we calculate instead a rough approximation to log |y log(x)|

   RR t1;
   long k;

   if (x > 0.5 && x < 1.5) { 
      xcopy(t1, x - 1);
      k = Lg2(t1);
   }
   else {
      k = NumBits(Lg2(x)); 
   }

   k += Lg2(y);

   if (k > 1000) Error("RR: overflow");

   if (k < 0) k = 0;

   
   RR::SetPrecision(p + k + 10);

   xcopy(t1, exp(y*log(x)));

   RR::SetPrecision(p);

   xcopy(res, t1);
}


void ReallyComputePi(RR& res)
{
   long p = RR::precision();
   RR::SetPrecision(p + NumBits(p) + 10);


   RR sum1;

   RR s, s1, t, t1;

   s = 0;
   t = 0.5;
   t1 = 0.5;

   long i;

   for (i = 3; ; i+=2) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t1, t1, -0.25);
      div(t, t1, i);
   }

   xcopy(sum1, s);


   RR g;

   inv(g, to_RR(3)); // g = 1/3

   s = 0;
   xcopy(t, g);
   xcopy(t1, g);

   sqr(g, g);
   negate(g, g); // g = -1/9

   for (i = 3; ; i+=2) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t1, t1, g);
      div(t, t1, i);
   }

   add(s, s, sum1);
   mul(s, s, 4);

   RR::SetPrecision(p);
   xcopy(res, s);
}

void ComputePi(RR& res)
{
   static long prec = 0;
   static RR pi;

   long p = RR::precision();

   if (prec <= p + 10) {
      prec = p + 20;
      RR::SetPrecision(prec);
      ReallyComputePi(pi);
      RR::SetPrecision(p);
   }

   xcopy(res, pi);
}



void sin(RR& res, const RR& x)
{
   if (x == 0) {
      res = 0;
      return;
   }

   if (Lg2(x) > 1000) 
      Error("sin: sorry...argument too large in absolute value");

   long p = RR::precision();

   RR pi, t1, f;
   RR n;

   // we want to make x^2 < 3, so that the series for sin(x)
   // converges nicely, without any nasty cancellations in the
   // first terms of the series.

   if (x*x < 3) {
      RR::SetPrecision(p + NumBits(p) + 10);
      xcopy(f, x);
   }
   else {

      // we want to write x/pi = n + f, |f| < 1/2....
      // but we have to do *this* very carefully, so that f is computed
      // to precision > p.  I know, this is sick!
   
      long p1;
   
      p1 = p + Lg2(x) + 20;
   
   
      for (;;) {
         RR::SetPrecision(p1);
         ComputePi(pi);
         xcopy(t1, x/pi);
         xcopy(n, floor(t1));
         xcopy(f, t1 - n);
         if (f > 0.5) {
            n++;
            xcopy(f, t1 - n);
         }

         if (f == 0 || p1 < p - Lg2(f) + Lg2(n) + 10) {
            // we don't have enough bits of f...increase p1 and continue

            p1 = p1 + max(20, p1/10);
         }
         else
            break;
      }

      RR::SetPrecision(p + NumBits(p) + 10);
      ComputePi(pi);

      xcopy(f, pi * f);
      
      if (n != 0 && n.exponent() == 0) {
         // n is odd, so we negate f, which negates sin(f)

         xcopy(f, -f);
      }

   }

   // Boy, that was painful, but now its over, and we can simply apply
   // the series for sin(f)

   RR t2, s, s1, t;
   long i;

   s = 0;
   xcopy(t, f);

   for (i = 3; ; i=i+2) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t, t, f);
      mul(t, t, f);
      div(t, t, i-1);
      div(t, t, i);
      negate(t, t);
   }

   RR::SetPrecision(p);
   
   xcopy(res, s);

}

void cos(RR& res, const RR& x)
{
   if (x == 0) {
      res = 1;
      return;
   }

   if (Lg2(x) > 1000) 
      Error("cos: sorry...argument too large in absolute value");

   long p = RR::precision();

   RR pi, t1, f;
   RR n;

   // we want to write x/pi = (n+1/2) + f, |f| < 1/2....
   // but we have to do *this* very carefully, so that f is computed
   // to precision > p.  I know, this is sick!

   long p1;

   p1 = p + Lg2(x) + 20;


   for (;;) {
      RR::SetPrecision(p1);
      ComputePi(pi);
      xcopy(t1, x/pi);
      xcopy(n, floor(t1));
      xcopy(f, t1 - (n + 0.5));

      if (f == 0 || p1 < p - Lg2(f) + Lg2(n) + 10) {
         // we don't have enough bits of f...increase p1 and continue

         p1 = p1 + max(20, p1/10);
      }
      else
         break;
   }

   RR::SetPrecision(p + NumBits(p) + 10);
   ComputePi(pi);

   xcopy(f, pi * f);
   
   if (n == 0 || n.exponent() != 0) {
      // n is even, so we negate f, which negates sin(f)

      xcopy(f, -f);
   }

   // Boy, that was painful, but now its over, and we can simply apply
   // the series for sin(f)

   RR t2, s, s1, t;
   long i;

   s = 0;
   xcopy(t, f);

   for (i = 3; ; i=i+2) {
      add(s1, s, t);
      if (s == s1) break;
      xcopy(s, s1);
      mul(t, t, f);
      mul(t, t, f);
      div(t, t, i-1);
      div(t, t, i);
      negate(t, t);
   }

   RR::SetPrecision(p);
   
   xcopy(res, s);

}

