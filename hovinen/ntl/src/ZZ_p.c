

#include <NTL/ZZ_p.h>
#include <NTL/tools.h>
#include <NTL/FFT.h>


ZZ_pInfoT::ZZ_pInfoT(const ZZ& NewP)
{
   if (NewP <= 1) Error("ZZ_pContext: p must be > 1");

   ref_count = 1;
   p = NewP;
   size = p.size();
   ExtendedModulusSize = 2*size + (NTL_BITS_PER_LONG + NTL_NBITS - 1)/NTL_NBITS;
   initialized = 0;
   x = 0;
   u = 0;
   tbl = 0;
   tbl1 = 0;

   long i;
   for (i = 0; i < MAX_ZZ_p_TEMPS; i++)
      temps[i] = 0;

   temps_top = 0;
}



void ZZ_pInfoT::init()
{
   ZZ B, M, M1, MinusM;
   long n, i;
   long q, t;

   initialized = 1;

   sqr(B, p);

   LeftShift(B, B, NTL_FFTMaxRoot+NTL_FFTFudge);

   set(M);
   n = 0;
   while (M <= B) {
      UseFFTPrime(n);
      q = FFTPrime[n];
      n++;
      mul(M, M, q);
   }

   NumPrimes = n;
   MaxRoot = CalcMaxRoot(q);

   double fn = double(n);

   // I've re-calculated the error bounds...
   // Although the bounds in versions < 3.7 were incorrect,
   // on any currently existing machine, this bug is "benign".

   if (fn*(fn+32) > 0.5*NTL_FDOUBLE_PRECISION)
      Error("modulus too big");


   if (fn*(fn+32) > (0.5*NTL_FDOUBLE_PRECISION)/double(NTL_RADIX))
      QuickCRT = 0;
   else
      QuickCRT = 1;


   negate(MinusM, M);
   rem(MinusMModP, MinusM, p);

   CoeffModP.SetSize(n, p.size());

   if (!(x = (double *) malloc(n * (sizeof (double)))))
      Error("out of space");

   if (!(u = (long *) malloc(n * (sizeof (long)))))
      Error("out of space");

   for (i = 0; i < n; i++) {
      q = FFTPrime[i];

      div(M1, M, q);
      t = rem(M1, q);
      t = InvMod(t, q);
      mul(M1, M1, t);
      rem(CoeffModP[i], M1, p);
      x[i] = ((double) t)/((double) q);
      u[i] = t;
   }

   B = p;
   LeftShift(B, B, NTL_NBITS);
   mul(B, B, NumPrimes);
   mmsize = B.size();

#if (defined(NTL_SINGLE_MUL))
   tbl = (double **) malloc(NumPrimes * sizeof(double *));
   if (!tbl) Error("out of space");
   for (i = 0; i < NumPrimes; i++) {
      tbl[i] = (double *) malloc(p.size() * sizeof(double));
      if (!tbl[i]) Error ("out of space");
   }

   long t1;
   long j;

   for (i = 0; i < NumPrimes; i++) {
      q = FFTPrime[i];
      t = (((long)1) << NTL_NBITS) % q;
      t1 = 1;
      tbl[i][0] = (double) 1;
      for (j = 1; j < p.size(); j++) {
         t1 = MulMod(t1, t, q);
         tbl[i][j] = (double) t1;
      }
   }
#else
   tbl = 0;
#endif

#if (defined(NTL_TBL_REM) && !defined(NTL_SINGLE_MUL))
   tbl1 = (long **) malloc(NumPrimes * sizeof(long *));
   if (!tbl1) Error("out of space");
   for (i = 0; i < NumPrimes; i++) {
      tbl1[i] = (long *) malloc(p.size() * sizeof(long));
      if (!tbl1[i]) Error ("out of space");
   }

   long t1;
   long j;

   for (i = 0; i < NumPrimes; i++) {
      q = FFTPrime[i];
      t = (((long)1) << NTL_NBITS) % q;
      t1 = 1;
      tbl1[i][0] = 1;
      for (j = 1; j < p.size(); j++) {
         t1 = MulMod(t1, t, q);
         tbl1[i][j] = t1;
      }
   }
#else
   tbl1 = 0;
#endif
}

ZZ_pInfoT::~ZZ_pInfoT()
{
   long i;

   for (i = 0; i < MAX_ZZ_p_TEMPS; i++)
      if (temps[i]) delete temps[i];

   if (initialized) {
      free(x);
      free(u);
#if (defined(NTL_SINGLE_MUL))
      for (i = 0; i < NumPrimes; i++)
         free(tbl[i]);

      free(tbl);

#elif (defined(NTL_TBL_REM))
      for (i = 0; i < NumPrimes; i++)
         free(tbl1[i]);

      free(tbl1);
#endif
   }
}


ZZ_pInfoT *ZZ_pInfo = 0; 

typedef ZZ_pInfoT *ZZ_pInfoPtr;


static 
void CopyPointer(ZZ_pInfoPtr& dst, ZZ_pInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative ZZ_pContext ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      src->ref_count++;

      if (src->ref_count < 0) 
         Error("internal error: ZZ_pContext ref_count overflow");
   }

   dst = src;
}
   


void ZZ_p::init(const ZZ& p)
{
   ZZ_pContext c(p);
   c.restore();
}


ZZ_pContext::ZZ_pContext(const ZZ& p)
{
   ptr = new ZZ_pInfoT(p);
}

ZZ_pContext::ZZ_pContext(const ZZ_pContext& a)
{
   ptr = 0;
   CopyPointer(ptr, a.ptr);
}

ZZ_pContext& ZZ_pContext::operator=(const ZZ_pContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}


ZZ_pContext::~ZZ_pContext()
{
   CopyPointer(ptr, 0);
}

void ZZ_pContext::save()
{
   CopyPointer(ptr, ZZ_pInfo);
}

void ZZ_pContext::restore() const
{
   CopyPointer(ZZ_pInfo, ptr);
}



ZZ_pBak::~ZZ_pBak()
{
   if (MustRestore)
      CopyPointer(ZZ_pInfo, ptr);

   CopyPointer(ptr, 0);
}

void ZZ_pBak::save()
{
   MustRestore = 1;
   CopyPointer(ptr, ZZ_pInfo);
}


void ZZ_pBak::restore()
{
   MustRestore = 0;
   CopyPointer(ZZ_pInfo, ptr);
}


ZZ_pTemp::ZZ_pTemp()
{
   if (ZZ_pInfo->temps_top == MAX_ZZ_p_TEMPS)
      Error("ZZ_p temporary: out of temps");

   pos = ZZ_pInfo->temps_top;
   ZZ_pInfo->temps_top++;
}

ZZ_pTemp::~ZZ_pTemp()
{
   ZZ_pInfo->temps_top--;
}

ZZ_p& ZZ_pTemp::val() const
{
   if (!ZZ_pInfo->temps[pos]) 
      ZZ_pInfo->temps[pos] = new ZZ_p;

   return *(ZZ_pInfo->temps[pos]);
}




const ZZ_p& ZZ_p::zero()
{
   static ZZ_p z(ZZ_p_NoAlloc);
   return z;
}

ZZ_p::DivHandlerPtr ZZ_p::DivHandler = 0;

ZZ_p::ZZ_p()
{
   rep.SetSize(ModulusSize());
}
   

ZZ_p::ZZ_p(INIT_VAL_TYPE, const ZZ& a) 
{
   rep.SetSize(ModulusSize());
   conv(*this, a);
} 

ZZ_p::ZZ_p(INIT_VAL_TYPE, long a)
{
   rep.SetSize(ModulusSize());
   conv(*this, a);
}


void conv(ZZ_p& x, long a)
{
   if (a == 0)
      clear(x);
   else if (a == 1)
      set(x);
   else {
      static ZZ y;

      conv(y, a);
      conv(x, y);
   }
}

istream& operator>>(istream& s, ZZ_p& x)
{
   static ZZ y;

   s >> y;
   conv(x, y);

   return s;
}

void div(ZZ_p& x, const ZZ_p& a, const ZZ_p& b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val(); 

   inv(T, b);
   mul(x, a, T);
}

void inv(ZZ_p& x, const ZZ_p& a)
{
   if (InvModStatus(x.rep, a.rep, ZZ_p::modulus())) {
      if (IsZero(a.rep))
         Error("ZZ_p: division by zero");
      else if (ZZ_p::DivHandler)
         (*ZZ_p::DivHandler)(a);
      else
         Error("ZZ_p: division by non-invertible element");
   }
}

long operator==(const ZZ_p& a, long b)
{
   if (b == 0)
      return IsZero(a);

   if (b == 1)
      return IsOne(a);

   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   return a == T;
}



void add(ZZ_p& x, const ZZ_p& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   add(x, a, T);
}

void sub(ZZ_p& x, const ZZ_p& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   sub(x, a, T);
}

void sub(ZZ_p& x, long a, const ZZ_p& b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, a);
   sub(x, T, b);
}

void mul(ZZ_p& x, const ZZ_p& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   mul(x, a, T);
}

void div(ZZ_p& x, const ZZ_p& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   div(x, a, T);
}

void div(ZZ_p& x, long a, const ZZ_p& b)
{
   if (a == 1) {
      inv(x, b);
   }
   else {
      ZZ_pTemp TT; ZZ_p& T = TT.val();
      conv(T, a);
      div(x, T, b);
   }
}

