

#include <NTL/lzz_pE.h>

zz_pEInfoT::zz_pEInfoT(const zz_pX& NewP)
{
   ref_count = 1;

   build(p, NewP);

   power(cardinality, zz_p::modulus(), deg(NewP));
}




zz_pEInfoT *zz_pEInfo = 0; 



typedef zz_pEInfoT *zz_pEInfoPtr;


static 
void CopyPointer(zz_pEInfoPtr& dst, zz_pEInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative zz_pEContext ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      src->ref_count++;

      if (src->ref_count < 0) 
         Error("internal error: zz_pEContext ref_count overflow");
   }

   dst = src;
}
   



void zz_pE::init(const zz_pX& p)
{
   zz_pEContext c(p);
   c.restore();
}


zz_pEContext::zz_pEContext(const zz_pX& p)
{
   ptr = new zz_pEInfoT(p);
}

zz_pEContext::zz_pEContext(const zz_pEContext& a)
{
   ptr = 0;
   CopyPointer(ptr, a.ptr);
}

zz_pEContext& zz_pEContext::operator=(const zz_pEContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}


zz_pEContext::~zz_pEContext()
{
   CopyPointer(ptr, 0);
}

void zz_pEContext::save()
{
   CopyPointer(ptr, zz_pEInfo);
}

void zz_pEContext::restore() const
{
   CopyPointer(zz_pEInfo, ptr);
}



zz_pEBak::~zz_pEBak()
{
   if (MustRestore)
      CopyPointer(zz_pEInfo, ptr);

   CopyPointer(ptr, 0);
}

void zz_pEBak::save()
{
   MustRestore = 1;
   CopyPointer(ptr, zz_pEInfo);
}



void zz_pEBak::restore()
{
   MustRestore = 0;
   CopyPointer(zz_pEInfo, ptr);
}



const zz_pE& zz_pE::zero()
{
   static zz_pE z(zz_pE_NoAlloc);
   return z;
}


zz_pE::zz_pE()
{
   rep.rep.SetMaxLength(zz_pE::degree());
}



istream& operator>>(istream& s, zz_pE& x)
{
   zz_pX y;

   s >> y;
   conv(x, y);

   return s;
}

void div(zz_pE& x, const zz_pE& a, const zz_pE& b)
{
   zz_pE t;

   inv(t, b);
   mul(x, a, t);
}

void div(zz_pE& x, const zz_pE& a, long b)
{
   NTL_zz_pRegister(B);
   B = b;
   inv(B, B);
   mul(x, a, B);
}

void div(zz_pE& x, const zz_pE& a, const zz_p& b)
{
   NTL_zz_pRegister(B);
   B = b;
   inv(B, B);
   mul(x, a, B);
}

void div(zz_pE& x, long a, const zz_pE& b)
{
   zz_pE t;
   inv(t, b);
   mul(x, a, t);
}

void div(zz_pE& x, const zz_p& a, const zz_pE& b)
{
   zz_pE t;
   inv(t, b);
   mul(x, a, t);
}



void inv(zz_pE& x, const zz_pE& a)
{
   InvMod(x.rep, a.rep, zz_pE::modulus());
}

