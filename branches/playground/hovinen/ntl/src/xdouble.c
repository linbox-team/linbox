
#include <NTL/xdouble.h>

#include <NTL/IsFinite.h>

#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>



long xdouble::oprec = 10;

void xdouble::SetOutputPrecision(long p)
{
   if (p < 1) Error("xdouble: output precision too small");

   oprec = p;
}

void xdouble::normalize() 
{
   if (x == 0) 
      e = 0;
   else if (x > 0) {
      while (x < NTL_XD_HBOUND_INV) { x *= NTL_XD_BOUND; e--; }
      while (x > NTL_XD_HBOUND) { x *= NTL_XD_BOUND_INV; e++; }
   }
   else {
      while (x > -NTL_XD_HBOUND_INV) { x *= NTL_XD_BOUND; e--; }
      while (x < -NTL_XD_HBOUND) { x *= NTL_XD_BOUND_INV; e++; }
   }

   if (e >= (1L << (NTL_BITS_PER_LONG-4)))
      Error("xdouble: overflow");

   if (e <= -(1L << (NTL_BITS_PER_LONG-4)))
      Error("xdouble: underflow");
}
   


xdouble to_xdouble(double a)
{
   if (a == 0 || a == 1 || (a > 0 && a >= NTL_XD_HBOUND_INV && a <= NTL_XD_HBOUND)
       || (a < 0 && a <= -NTL_XD_HBOUND_INV && a >= -NTL_XD_HBOUND)) {
      
      return xdouble(a, 0); 

   }

   if (!IsFinite(&a))
      Error("double to xdouble conversion: non finite value");

   xdouble z = xdouble(a, 0);
   z.normalize();
   return z;
}


void conv(double& xx, const xdouble& a)
{
   double x;
   long e;

   x = a.x;
   e = a.e;

   while (e > 0) { x *= NTL_XD_BOUND; e--; }
   while (e < 0) { x *= NTL_XD_BOUND_INV; e++; }

   xx = x;
}




xdouble operator+(const xdouble& a, const xdouble& b)
{
   xdouble z;

   if (a.x == 0) 
      return b;

   if (b.x == 0)
     return a;
      

   if (a.e == b.e) {
      z.x = a.x + b.x;
      z.e = a.e;
      z.normalize();
      return z;
   }
   else if (a.e > b.e) {
      if (a.e > b.e+1)
         return a;

      z.x = a.x + b.x*NTL_XD_BOUND_INV;
      z.e = a.e;
      z.normalize();
      return z;
   }
   else {
      if (b.e > a.e+1)
         return b;

      z.x = a.x*NTL_XD_BOUND_INV + b.x;
      z.e = b.e;
      z.normalize();
      return z;
   }
}


xdouble operator-(const xdouble& a, const xdouble& b)
{
   xdouble z;

   if (a.x == 0)
      return -b;

   if (b.x == 0)
      return a;

   if (a.e == b.e) {
      z.x = a.x - b.x;
      z.e = a.e;
      z.normalize();
      return z;
   }
   else if (a.e > b.e) {
      if (a.e > b.e+1)
         return a;

      z.x = a.x - b.x*NTL_XD_BOUND_INV;
      z.e = a.e;
      z.normalize();
      return z;
   }
   else {
      if (b.e > a.e+1)
         return -b;

      z.x = a.x*NTL_XD_BOUND_INV - b.x;
      z.e = b.e;
      z.normalize();
      return z;
   }
}

xdouble operator-(const xdouble& a)
{
   xdouble z;
   z.x = -a.x;
   z.e = a.e;
   return z;
}

xdouble operator*(const xdouble& a, const xdouble& b)
{
   xdouble z;

   z.e = a.e + b.e;
   z.x = a.x * b.x;
   z.normalize();
   return z;
}

xdouble operator/(const xdouble& a, const xdouble& b)
{
   xdouble z;

   if (b.x == 0) Error("zdouble division by 0");

   z.e = a.e - b.e;
   z.x = a.x / b.x;
   z.normalize();
   return z;
}



long compare(const xdouble& a, const xdouble& b)
{
   xdouble z = a - b;

   if (z.x < 0)
      return -1;
   else if (z.x == 0)
      return 0;
   else
      return 1;
}

long sign(const xdouble& z)
{
   if (z.x < 0)
      return -1;
   else if (z.x == 0)
      return 0;
   else
      return 1;
}
   


xdouble trunc(const xdouble& a)
{
   if (a.x >= 0)
      return floor(a);
   else
      return ceil(a);
}


xdouble floor(const xdouble& a)
{
   xdouble z;

   if (a.e == 0) {
      z.x = floor(a.x);
      z.e = 0;
      z.normalize();
      return z;
   }
   else if (a.e > 0) {
      return a;
   }
   else {
      if (a.x < 0)
         return to_xdouble(-1);
      else
         return to_xdouble(0);
   }
}

xdouble ceil(const xdouble& a)
{
   xdouble z;

   if (a.e == 0) {
      z.x = ceil(a.x);
      z.e = 0;
      z.normalize();
      return z;
   }
   else if (a.e > 0) {
      return a;
   }
   else {
      if (a.x < 0)
         return to_xdouble(0);
      else
         return to_xdouble(1);
   }
}

xdouble to_xdouble(const ZZ& a)
{
   long *n = a.rep;
   xdouble res;
   long i;
   static xdouble fradix = to_xdouble(NTL_RADIX);

   if (!n) return to_xdouble(0);

   if ((i = n[0]) < 0)
      i = -i;
   res = to_xdouble(n[i--]);
   for (; i; i--)
      res = res * fradix + n[i];
   if (n[0] < 0)
      res = -res;

   return res;
}



#define MustAlloc(c, len)  (!(c) || ((c)[-1] >> 1) < (len))



void zxdoubtoz(const xdouble& aa, verylong *xx)
{
   verylong x;
   long neg, i, t, sz;
   static xdouble fradix_inv = to_xdouble(NTL_FRADIX_INV);
   static xdouble fradix = to_xdouble(NTL_RADIX);


   xdouble a;

   a = floor(aa);

   ForceToMem(&a.x);

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
   a = a*fradix_inv;

   while (a >= 1) {
      a = a*fradix_inv;
      sz++;
   }

   x = *xx;
   if (MustAlloc(x, sz)) {
      zsetlength(&x, sz);
      *xx = x;
   }


   for (i = sz; i > 0; i--) {
      a = a*fradix;
      conv(t, a);
      x[i] = t;
      a = a - t;
   }

   x[0] = (neg ? -sz : sz);
}

xdouble fabs(const xdouble& a)
{
   xdouble z;

   z.e = a.e;
   z.x = fabs(a.x);
   return z;
}

xdouble sqrt(const xdouble& a)
{
   if (a == 0)
      return to_xdouble(0);

   if (a < 0)
      Error("xdouble: sqrt of negative number");

   xdouble t;

   if (a.e & 1) {
      t.e = (a.e - 1)/2;
      t.x = sqrt(a.x * NTL_XD_BOUND);
   }
   else {
      t.e = a.e/2;
      t.x = sqrt(a.x);
   }

   t.normalize();

   return t;
}
      

void power(xdouble& z, const xdouble& a, const ZZ& e)
{
   xdouble b, res;

   b = a;

   res = 1;
   long n = NumBits(e);
   long i;

   for (i = n-1; i >= 0; i--) {
      res = res * res;
      if (bit(e, i))
         res = res * b;
   }

   if (sign(e) < 0) 
      z = 1/res;
   else
      z = res;
}




void power(xdouble& z, const xdouble& a, long e)
{
   static ZZ E;
   E = e;
   power(z, a, E);
}
   

   

ostream& operator<<(ostream& s, const xdouble& a)
{
   if (a == 0) {
      s << "0";
      return s;
   }

   xdouble b;
   long neg;

   if (a < 0) {
      b = -a;
      neg = 1;
   }
   else {
      b = a;
      neg = 0;
   }

   long k;

   k = long((xdouble::oprec*3.321928095-log(b.x)/log(2.0)-b.e)/3.321928095);


   xdouble c;

   power(c, to_xdouble(10), k);

   b = b*c;

   power(c, to_xdouble(10), xdouble::oprec);

   while (b < c) {
      b = b * 10;
      k++;
   }

   while (b >= c) {
      b = b/10;
      k--;
   }

   b = b + 0.5;
   k = -k;

   ZZ B;
   conv(B, b);

   char *bp = new char[xdouble::oprec+10];

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

   delete [] bp;
   return s;
}

istream& operator>>(istream& s, xdouble& x)
{
   long c;
   long sign;
   ZZ a, b;

   if (!s) Error("bad xdouble input");

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

      if (c < '0' || c > '9') Error("bad xdouble input");

      e = 0;
      while (c >= '0' && c <= '9') {
         mul(e, e, 10);
         add(e, e, c-'0');
         s.get();
         c = s.peek();
      }
   }

   if (!got1 && !got2 && !got_e) Error("bad xdouble input");

   xdouble t1, t2, v;

   if (got1 || got2) {
      conv(t1, a);
      conv(t2, b);
      v = t1/t2;
   }
   else
      v = 1;

   if (sign < 0)
      v = -v;

   if (got_e) {
      if (e_sign < 0) negate(e, e);
      power(t1, to_xdouble(10), e);
      v = v * t1;
   }

   x = v;
   return s;
}
      

void power2(xdouble& z, long e)
{
   power(z, to_xdouble(2), e);
}


void MulAdd(xdouble& z, const xdouble& a, const xdouble& b, const xdouble& c)
// z = a + b*c
{
   double x;
   long e;

   e = b.e + c.e;
   x = b.x * c.x;

   if (x == 0) { 
      z = a;
      return;
   }

   if (a.x == 0) {
      z.e = e;
      z.x = x;
      z.normalize();
      return;
   }
      

   if (a.e == e) {
      z.x = a.x + x;
      z.e = e;
      z.normalize();
      return;
   }
   else if (a.e > e) {
      if (a.e > e+1) {
         z = a;
         return;
      }

      z.x = a.x + x*NTL_XD_BOUND_INV;
      z.e = a.e;
      z.normalize();
      return;
   }
   else {
      if (e > a.e+1) {
         z.x = x;
         z.e = e;
         z.normalize();
         return;
      }

      z.x = a.x*NTL_XD_BOUND_INV + x;
      z.e = e;
      z.normalize();
      return;
   }
}

void MulSub(xdouble& z, const xdouble& a, const xdouble& b, const xdouble& c)
// z = a - b*c
{
   double x;
   long e;

   e = b.e + c.e;
   x = b.x * c.x;

   if (x == 0) { 
      z = a;
      return;
   }

   if (a.x == 0) {
      z.e = e;
      z.x = -x;
      z.normalize();
      return;
   }
      

   if (a.e == e) {
      z.x = a.x - x;
      z.e = e;
      z.normalize();
      return;
   }
   else if (a.e > e) {
      if (a.e > e+1) {
         z = a;
         return;
      }

      z.x = a.x - x*NTL_XD_BOUND_INV;
      z.e = a.e;
      z.normalize();
      return;
   }
   else {
      if (e > a.e+1) {
         z.x = -x;
         z.e = e;
         z.normalize();
         return;
      }

      z.x = a.x*NTL_XD_BOUND_INV - x;
      z.e = e;
      z.normalize();
      return;
   }
}

double log(const xdouble& a)
{
   static double LogBound = log(NTL_XD_BOUND);
   if (a.x <= 0) {
      Error("log(xdouble): argument must be positive");
   }

   return log(a.x) + a.e*LogBound;
}

xdouble xexp(double x)
{
   const double LogBound = log(NTL_XD_BOUND);

   double y = x/LogBound;
   double iy = floor(y+0.5);

   if (iy >= (1L << (NTL_BITS_PER_LONG-4)))
      Error("xdouble: overflow");

   if (iy <= -(1L << (NTL_BITS_PER_LONG-4)))
      Error("xdouble: underflow");


   double fy = y - iy;

   xdouble res;
   res.e = long(iy);
   res.x = exp(fy*LogBound);
   res.normalize();
   return res;
}

xdouble to_xdouble(const char *s)
{
   long c;
   long sign;
   ZZ a, b;
   long i=0;

   if (!s) Error("bad xdouble input");

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

      if (c < '0' || c > '9') Error("bad xdouble input");

      e = 0;
      while (c >= '0' && c <= '9') {
         mul(e, e, 10);
         add(e, e, c-'0');
         i++;
         c = s[i];
      }
   }

   if (!got1 && !got2 && !got_e) Error("bad xdouble input");

   xdouble t1, t2, v;

   if (got1 || got2) {
      conv(t1, a);
      conv(t2, b);
      v = t1/t2;
   }
   else
      v = 1;

   if (sign < 0)
      v = -v;

   if (got_e) {
      if (e_sign < 0) negate(e, e);
      power(t1, to_xdouble(10), e);
      v = v * t1;
   }

   return v;
}

