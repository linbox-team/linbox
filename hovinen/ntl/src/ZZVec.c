
#include <NTL/ZZVec.h>
#include <NTL/tools.h>



void ZZVec::SetSize(long n, long d)
{
   if (n < 0 || d <= 0) Error("bad args to ZZVec::SetSize()");

   if (v)
      Error("ZZVec initialized more than once");

   if (d >= (1L << (NTL_BITS_PER_LONG-4))/NTL_NBITS)
      Error("size too big in ZZVec::SetSize");

   if (n == 0) return;

   long size = (d+3);
   long AllocAmt = MaxAllocBlock / size;
   if (AllocAmt == 0) AllocAmt = 1;

   v = (ZZ*) malloc(n * (sizeof (ZZ)));
   if (!v) Error("out of memory in ZZVec::SetSize()");

   long i = 0;
   long m;
   long *p, *q;
   long j;

   while (i < n) {
      m = min((n-i), AllocAmt);
      p = (long *) malloc(m*size*(sizeof (long)));
      if (!p) Error("out of memory in ZZVec::SetSize()");
      for (j = 0, q = p+1; j < m; j++, i++, q += size) {
         q[-1] = ((d+1) << 1) | 1;;
         q[0] = 1;
         q[1] = 0;
         v[i].rep = q; 
      }
   }

   len = n;
   bsize = d;
}

void ZZVec::kill()
{
   long n = len;
   long d = bsize;
   if (n == 0) return;

   long size = (d+3);
   long AllocAmt = MaxAllocBlock / size;
   if (AllocAmt == 0) AllocAmt = 1;

   long i = 0;
   long m;

   while (i < n) {
      m = min((n-i), AllocAmt);
      free(v[i].rep-1);
      i += m;
   }

   free(v);

   v = 0; len = 0; bsize = 0;
}


ZZVec& ZZVec::operator=(const ZZVec& a) 
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
   
ZZVec::ZZVec(const ZZVec& a)
{
   v = 0; len = 0; bsize = 0;

   SetSize(a.len, a.bsize);

   long i;
   for (i = 0; i < a.len; i++)
      v[i] = (a.v)[i];
}

void swap(ZZVec& x, ZZVec& y)
{
   ZZ* t1;
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

