
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <math.h>

xdouble xpow(double a, double b)
{
   return xexp(b*log(a));
}

double log2_pkt(long kk, long tt)
{
   double k = kk;
   double t = tt;

   double e0;
   e0 = -t*log(4.0);
   // This bound follows from Burthe's PhD thesis (1995).
   // See Bach & Shallit's book.

   double e1;
   e1 = log(xpow(k, 2) * xpow(4, 2-sqrt(k)));

   double e2;
   if ((t == 2 && k >= 88) || (3 <= t && t <= k/9 && k >= 21))
      e2 = log(xpow(k, 1.5) * xpow(2, t) * xpow(t, -0.5) * xpow(4, 2 - sqrt(t*k)));
   else 
      e2 = 0;

   double e3;

   if (k/9 <= t && t <= k/4 && k >= 21) {
      e3 = log( (7.0/20.0)*k*xpow(2, -5*t) +
                (1.0/7.0)*xpow(k, 15.0/4.0)*xpow(2, -k/2-2*t) +
                12*k*xpow(2, -k/4-3*t) );
   }
   else
      e3 = 0;

   double e4;
   if (t >= k/4 && k >= 21)
      e4 = log( (1.0/7.0)*xpow(k, 15.0/4.0)*xpow(2, -k/2-2*t) );
   else
      e4 = 0;

   double emin;

   emin = e0;
   if (e1 < emin) emin = e1;
   if (e2 < emin) emin = e2;
   if (e3 < emin) emin = e3;
   if (e4 < emin) emin = e4;

   return emin/log(2.0);
}


void GenPrime(ZZ& n, long k, long err)
{
   if (k <= 1) Error("GenPrime: bad length");

   if (k == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }

   long t;

   for (t = 1; ; t++) {
      if (log2_pkt(k, t) <= -err)
         break;
   }

   RandomPrime(n, k, t);
}


long GenPrime_long(long k, long err)
{
   if (k <= 1) Error("GenPrime: bad length");

   if (k == 2) {
      if (RandomBnd(2))
         return 3;
      else
         return 2;
   }

   long t;

   for (t = 1; ; t++) {
      if (log2_pkt(k, t) <= -err)
         break;
   }

   return RandomPrime_long(k, t);
}

