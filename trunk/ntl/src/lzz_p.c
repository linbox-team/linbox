
#include <NTL/lzz_p.h>
#include <NTL/FFT.h>

zz_pInfoT::zz_pInfoT(long NewP, long maxroot)
{
   ref_count = 1;

   if (NewP <= 1) Error("zz_pContext: p must be > 1");
   if (NumBits(NewP) > NTL_NBITS) Error("zz_pContext: modulus too big");

   ZZ P, B, M, M1, MinusM;
   long n, i;
   long q, t;

   p = NewP;

   pinv = 1/double(p);

   index = -1;

   ::conv(P, p);

   ::sqr(B, P);
   LeftShift(B, B, maxroot+NTL_FFTFudge);

   set(M);
   n = 0;
   while (M <= B) {
      UseFFTPrime(n);
      q = FFTPrime[n];
      n++;
      ::mul(M, M, q);
   }

   if (n > 4) Error("zz_pInit: too many primes");

   NumPrimes = n;
   PrimeCnt = n;
   MaxRoot = CalcMaxRoot(q);

   if (maxroot < MaxRoot)
      MaxRoot = maxroot;

   ::negate(MinusM, M);
   MinusMModP = rem(MinusM, p);

   if (!(CoeffModP = (long *) malloc(n * (sizeof (long)))))
      Error("out of space");

   if (!(x = (double *) malloc(n * (sizeof (double)))))
      Error("out of space");

   if (!(u = (long *) malloc(n * (sizeof (long)))))
      Error("out of space");

   for (i = 0; i < n; i++) {
      q = FFTPrime[i];

      ::div(M1, M, q);
      t = rem(M1, q);
      t = InvMod(t, q);
      ::mul(M1, M1, t);
      CoeffModP[i] = rem(M1, p);
      x[i] = ((double) t)/((double) q);
      u[i] = t;
   }
}

zz_pInfoT::zz_pInfoT(long Index)
{
   ref_count = 1;

   index = Index;

   if (index < 0)
      Error("bad FFT prime index");

   // allows non-consecutive indices...I'm not sure why
   while (NumFFTPrimes < index)
      UseFFTPrime(NumFFTPrimes);

   UseFFTPrime(index);

   p = FFTPrime[index];
   pinv = FFTPrimeInv[index];

   NumPrimes = 1;
   PrimeCnt = 0;

   MaxRoot = CalcMaxRoot(p);
}




zz_pInfoT::~zz_pInfoT()
{
   if (index < 0) {
      free(CoeffModP);
      free(x);
      free(u);
   }
}


zz_pInfoT *zz_pInfo = 0;


typedef zz_pInfoT *zz_pInfoPtr;

#if ((defined (_THREAD_SAFE)) || (defined (_REENTRANT))) \
      && (!defined (COARSE_LOCKS))

pthread_rwlock_t zz_p_lock;

#endif


static 
void CopyPointer(zz_pInfoPtr& dst, zz_pInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative zz_pContext ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      src->ref_count++;

      if (src->ref_count < 0) 
         Error("internal error: zz_pContext ref_count overflow");
   }

   dst = src;
}
   

void zz_p::init(long p, long maxroot)
{
   zz_pContext c(p, maxroot);
   c.restore();

}

void zz_p::FFTInit(long index)
{
   zz_pContext c(INIT_FFT, index);
   c.restore();
}

zz_pContext::zz_pContext(long p, long maxroot)
{
   ptr = new zz_pInfoT(p, maxroot);
}

zz_pContext::zz_pContext(INIT_FFT_TYPE, long index)
{
   ptr = new zz_pInfoT(index);
}

zz_pContext::zz_pContext(const zz_pContext& a)
{
   ptr = 0;
   CopyPointer(ptr, a.ptr);
}


zz_pContext& zz_pContext::operator=(const zz_pContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}


zz_pContext::~zz_pContext()
{
   CopyPointer(ptr, 0);
}

void zz_pContext::save()
{
   CopyPointer(ptr, zz_pInfo);
}

void zz_pContext::restore() const
{
#if (defined (_THREAD_SAFE)) || (defined (_REENTRANT))
   pthread_rwlock_wrlock (&zz_p_lock);
#endif

   CopyPointer(zz_pInfo, ptr);
   
#if (defined (_THREAD_SAFE)) || (defined (_REENTRANT))
   pthread_rwlock_unlock (&zz_p_lock);
#endif
}



zz_pBak::~zz_pBak()
{
   if (MustRestore)
      CopyPointer(zz_pInfo, ptr);

   CopyPointer(ptr, 0);
}

void zz_pBak::save()
{
   MustRestore = 1;
   CopyPointer(ptr, zz_pInfo);
}


void zz_pBak::restore()
{
   MustRestore = 0;

#if (defined (_THREAD_SAFE)) || (defined (_REENTRANT))
   pthread_rwlock_wrlock (&zz_p_lock);
#endif

   CopyPointer(zz_pInfo, ptr);
   
#if (defined (_THREAD_SAFE)) || (defined (_REENTRANT))
   pthread_rwlock_unlock (&zz_p_lock);
#endif
}




inline long reduce(long a, long p)
{
   if (a >= 0 && a < p)
      return a;
   else {
      a = a % p;
      if (a < 0) a += p;
      return a;
   }
}


zz_p zz_pInfoT::to_zz_p(long a) const
{
   return zz_p(reduce(a, p), INIT_LOOP_HOLE);
}

void zz_pInfoT::conv(zz_p& x, long a) const
{
   x.rep = reduce(a, p);
}

zz_p zz_pInfoT::to_zz_p(const ZZ& a) const
{
   return zz_p(rem(a, p), INIT_LOOP_HOLE);
}

void zz_pInfoT::conv(zz_p& x, const ZZ& a) const
{
   x.rep = rem(a, p);
}


istream& operator>>(istream& s, zz_p& x)
{
   _BUFFER ZZ y;
   s >> y;
   conv(x, y);

   return s;
}

ostream& operator<<(ostream& s, zz_p a)
{
   _BUFFER ZZ y;
   y = rep(a);
   s << y;

   return s;
}
