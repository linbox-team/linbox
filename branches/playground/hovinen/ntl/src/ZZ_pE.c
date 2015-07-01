

#include <NTL/ZZ_pE.h>

ZZ_pEInfoT::ZZ_pEInfoT(const ZZ_pX& NewP)
{
   ref_count = 1;

   build(p, NewP);

   power(cardinality, ZZ_p::modulus(), deg(NewP));
}




ZZ_pEInfoT *ZZ_pEInfo = 0; 



typedef ZZ_pEInfoT *ZZ_pEInfoPtr;


static 
void CopyPointer(ZZ_pEInfoPtr& dst, ZZ_pEInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative ZZ_pEContext ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      src->ref_count++;

      if (src->ref_count < 0) 
         Error("internal error: ZZ_pEContext ref_count overflow");
   }

   dst = src;
}
   



void ZZ_pE::init(const ZZ_pX& p)
{
   ZZ_pEContext c(p);
   c.restore();
}


ZZ_pEContext::ZZ_pEContext(const ZZ_pX& p)
{
   ptr = new ZZ_pEInfoT(p);
}

ZZ_pEContext::ZZ_pEContext(const ZZ_pEContext& a)
{
   ptr = 0;
   CopyPointer(ptr, a.ptr);
}

ZZ_pEContext& ZZ_pEContext::operator=(const ZZ_pEContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}


ZZ_pEContext::~ZZ_pEContext()
{
   CopyPointer(ptr, 0);
}

void ZZ_pEContext::save()
{
   CopyPointer(ptr, ZZ_pEInfo);
}

void ZZ_pEContext::restore() const
{
   CopyPointer(ZZ_pEInfo, ptr);
}



ZZ_pEBak::~ZZ_pEBak()
{
   if (MustRestore)
      CopyPointer(ZZ_pEInfo, ptr);

   CopyPointer(ptr, 0);
}

void ZZ_pEBak::save()
{
   MustRestore = 1;
   CopyPointer(ptr, ZZ_pEInfo);
}



void ZZ_pEBak::restore()
{
   MustRestore = 0;
   CopyPointer(ZZ_pEInfo, ptr);
}



const ZZ_pE& ZZ_pE::zero()
{
   static ZZ_pE z(ZZ_pE_NoAlloc);
   return z;
}


ZZ_pE::ZZ_pE()
{
   rep.rep.SetMaxLength(ZZ_pE::degree());
}



istream& operator>>(istream& s, ZZ_pE& x)
{
   ZZ_pX y;

   s >> y;
   conv(x, y);

   return s;
}

void div(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
{
   ZZ_pE t;

   inv(t, b);
   mul(x, a, t);
}

void div(ZZ_pE& x, const ZZ_pE& a, long b)
{
   NTL_ZZ_pRegister(B);
   B = b;
   inv(B, B);
   mul(x, a, B);
}

void div(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
{
   NTL_ZZ_pRegister(B);
   B = b;
   inv(B, B);
   mul(x, a, B);
}

void div(ZZ_pE& x, long a, const ZZ_pE& b)
{
   ZZ_pE t;
   inv(t, b);
   mul(x, a, t);
}

void div(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
{
   ZZ_pE t;
   inv(t, b);
   mul(x, a, t);
}



void inv(ZZ_pE& x, const ZZ_pE& a)
{
   InvMod(x.rep, a.rep, ZZ_pE::modulus());
}

