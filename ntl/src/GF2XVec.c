
#include <NTL/GF2XVec.h>
#include <NTL/tools.h>



void GF2XVec::SetSize(long n, long d)
{
   if (n < 0 || d <= 0) Error("bad args to GF2XVec::SetSize()");
   if (v)
      Error("GF2XVec initialized more than once");

   if (d >= (1L << (NTL_BITS_PER_LONG-4))/NTL_BITS_PER_LONG)
      Error("size too big in GF2XVec::SetSize");

   if (n == 0) return;

   long size = (d+2);
   long AllocAmt = MaxAllocBlock / size;
   if (AllocAmt == 0) AllocAmt = 1;

   v = (GF2X*) malloc(n * (sizeof (GF2X)));
   if (!v) Error("out of memory in GF2XVec::SetSize()");

   long i = 0;
   long m;
   u_long *p, *q;
   long j;

   while (i < n) {
      m = min((n-i), AllocAmt);
      p = (u_long *) malloc(m*size*(sizeof (u_long)));
      if (!p) Error("out of memory in GF2XVec::SetSize()");
      for (j = 0, q = p+2; j < m; j++, i++, q += size) {
         q[-2] = (d << 1) | 1;
         q[-1] = 0;
         v[i].xrep.rep = q; 
      }
   }

   len = n;
   bsize = d;
}

void GF2XVec::kill()
{
   long n = len;
   long d = bsize;
   if (n == 0) return;

   long size = (d+2);
   long AllocAmt = MaxAllocBlock / size;
   if (AllocAmt == 0) AllocAmt = 1;

   long i = 0;
   long m;

   while (i < n) {
      m = min((n-i), AllocAmt);
      free(v[i].xrep.rep-2);
      i += m;
   }

   free(v);

   v = 0; len = 0; bsize = 0;
}


GF2XVec& GF2XVec::operator=(const GF2XVec& a) 
{
   if (this == &a)
      return *this;

   kill();
   SetSize(a.len, a.bsize);

   long i;
   for (i = 0; i < a.len; i++)
      v[i] = (a.v)[i];

   return *this;
}
   
GF2XVec::GF2XVec(const GF2XVec& a)
{
   v = 0; len = 0; bsize = 0;

   SetSize(a.len, a.bsize);

   long i;
   for (i = 0; i < a.len; i++)
      v[i] = (a.v)[i];
}

void swap(GF2XVec& x, GF2XVec& y)
{
   GF2X* t1;
   long t2;

   t1 = x.v;
   x.v = y.v;
   y.v = t1;

   t2 = x.len;
   x.len = y.len;
   y.len = t2;

   t2 = x.bsize;
   x.bsize = y.bsize;
   y.bsize = t2;
}

