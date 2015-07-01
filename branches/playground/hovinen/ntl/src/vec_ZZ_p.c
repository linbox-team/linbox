

#include <NTL/vec_ZZ_p.h>
#include <NTL/tools.h>



void BlockConstruct(ZZ_p* x, long n)
{
   if (n <= 0) return; 

   long d = ZZ_p::ModulusSize();

   long size = d  + 3;
   long AllocAmt = (MaxAllocBlock-1) / size;
   if (AllocAmt == 0) AllocAmt = 1;

   long i = 0;
   long m;
   long *p, *q;
   long j;

   while (i < n) {
      m = min((n-i), AllocAmt);
      p = (long *) malloc((m*size + 1)*(sizeof (long)));
      if (!p) Error("out of memory in BlockConstruct(ZZ_p*,long)");
      *p = m;
      for (j = 0, q = p+2; j < m; j++, i++, q += size) {
         q[-1] = ((d+1) << 1 ) | 1;
         q[0] = 1;
         q[1] = 0;
         x[i].rep.rep = q;
      }
   }
}

void BlockDestroy(ZZ_p* x, long n)
{
   if (n <= 0) return;

   long i = 0;
   long *p;
   long m;

   while (i < n) {
      p = x[i].rep.rep-2;
      m =  *p;
      free(p);
      i += m;
   }
}


NTL_vector_impl_plain(ZZ_p,vec_ZZ_p)

NTL_io_vector_impl(ZZ_p,vec_ZZ_p)

NTL_eq_vector_impl(ZZ_p,vec_ZZ_p)


void conv(vec_ZZ_p& x, const vec_ZZ& a)
{
   long i, n;

   n = a.length();
   x.SetLength(n);

   ZZ_p* xp = x.elts();
   const ZZ* ap = a.elts();

   for (i = 0; i < n; i++)
      conv(xp[i], ap[i]);
}

void conv(vec_ZZ& x, const vec_ZZ_p& a)
{
   long n = a.length();
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      x[i] = rep(a[i]);
}



void InnerProduct(ZZ_p& x, const vec_ZZ_p& a, const vec_ZZ_p& b)
{
   long n = min(a.length(), b.length());
   long i;
   static ZZ accum, t;

   clear(accum);
   for (i = 0; i < n; i++) {
      mul(t, rep(a[i]), rep(b[i]));
      add(accum, accum, t);
   }

   conv(x, accum);
}

void InnerProduct(ZZ_p& x, const vec_ZZ_p& a, const vec_ZZ_p& b,
                  long offset)
{
   long n = min(a.length(), b.length()+offset);
   long i;
   static ZZ accum, t;

   clear(accum);
   for (i = offset; i < n; i++) {
      mul(t, rep(a[i]), rep(b[i-offset]));
      add(accum, accum, t);
   }

   conv(x, accum);
}

void mul(vec_ZZ_p& x, const vec_ZZ_p& a, const ZZ_p& b_in)
{
   NTL_ZZ_pRegister(b);
   b = b_in;
   long n = a.length();
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      mul(x[i], a[i], b);
}

void mul(vec_ZZ_p& x, const vec_ZZ_p& a, long b_in)
{
   NTL_ZZ_pRegister(b);
   b = b_in;
   long n = a.length();
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      mul(x[i], a[i], b);
}


void add(vec_ZZ_p& x, const vec_ZZ_p& a, const vec_ZZ_p& b)
{
   long n = a.length();
   if (b.length() != n) Error("vector add: dimension mismatch");

   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      add(x[i], a[i], b[i]);
}

void sub(vec_ZZ_p& x, const vec_ZZ_p& a, const vec_ZZ_p& b)
{
   long n = a.length();
   if (b.length() != n) Error("vector sub: dimension mismatch");
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      sub(x[i], a[i], b[i]);
}

void clear(vec_ZZ_p& x)
{
   long n = x.length();
   long i;
   for (i = 0; i < n; i++)
      clear(x[i]);
}

void negate(vec_ZZ_p& x, const vec_ZZ_p& a)
{
   long n = a.length();
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      negate(x[i], a[i]);
}

long IsZero(const vec_ZZ_p& a)
{
   long n = a.length();
   long i;

   for (i = 0; i < n; i++)
      if (!IsZero(a[i]))
         return 0;

   return 1;
}

vec_ZZ_p operator+(const vec_ZZ_p& a, const vec_ZZ_p& b)
{
   vec_ZZ_p res;
   add(res, a, b);
   NTL_OPT_RETURN(vec_ZZ_p, res);
}

vec_ZZ_p operator-(const vec_ZZ_p& a, const vec_ZZ_p& b)
{
   vec_ZZ_p res;
   sub(res, a, b);
   NTL_OPT_RETURN(vec_ZZ_p, res);
}


vec_ZZ_p operator-(const vec_ZZ_p& a)
{
   vec_ZZ_p res;
   negate(res, a);
   NTL_OPT_RETURN(vec_ZZ_p, res);
}


ZZ_p operator*(const vec_ZZ_p& a, const vec_ZZ_p& b)
{
   ZZ_p res;
   InnerProduct(res, a, b);
   NTL_OPT_RETURN(ZZ_p, res);
}


void VectorCopy(vec_ZZ_p& x, const vec_ZZ_p& a, long n)
{
   if (n < 0) Error("VectorCopy: negative length");
   if (n >= (1L << (NTL_BITS_PER_LONG-4))) Error("overflow in VectorCopy");

   long m = min(n, a.length());

   x.SetLength(n);
   
   long i;

   for (i = 0; i < m; i++)
      x[i] = a[i];

   for (i = m; i < n; i++)
      clear(x[i]);
}


