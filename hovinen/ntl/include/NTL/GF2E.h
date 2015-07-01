

#ifndef NTL_GF2E__H
#define NTL_GF2E__H

#include <NTL/GF2X.h>



class GF2EInfoT {
private:
   GF2EInfoT();                       // disabled
   GF2EInfoT(const GF2EInfoT&);   // disabled
   void operator=(const GF2EInfoT&);  // disabled
public:
   long ref_count;

   GF2EInfoT(const GF2X& NewP);
   ~GF2EInfoT() { }

   GF2XModulus p;

   long KarCross;
   long ModCross;
   long DivCross;

   ZZ cardinality;
};

extern GF2EInfoT *GF2EInfo; // info for current modulus, initially null




class GF2EContext {
private:
GF2EInfoT *ptr;

public:
void save();
void restore() const;

GF2EContext() { ptr = 0; }
GF2EContext(const GF2X& p);

GF2EContext(const GF2EContext&); 


GF2EContext& operator=(const GF2EContext&); 


~GF2EContext();


};


class GF2EBak {
private:
long MustRestore;
GF2EInfoT *ptr;

GF2EBak(const GF2EBak&); // disabled
void operator=(const GF2EBak&); // disabled

public:
void save();
void restore();

GF2EBak() { MustRestore = 0; ptr = 0; }

~GF2EBak();


};



struct GF2E_NoAlloc_type { GF2E_NoAlloc_type() { } };
const GF2E_NoAlloc_type GF2E_NoAlloc = GF2E_NoAlloc_type();



class GF2E {


private:

GF2X rep;

// static data


public:


// ****** constructors and assignment

GF2E() { rep.xrep.SetMaxLength(GF2E::WordLength()); }

GF2E(GF2E& x, INIT_TRANS_TYPE) : rep(x.rep, INIT_TRANS) { }

GF2E(const GF2E& a)  
   { rep.xrep.SetMaxLength(GF2E::WordLength()); rep = a.rep; }

GF2E(GF2E_NoAlloc_type) { }  // allocates no space

~GF2E() { } 

GF2E& operator=(const GF2E& a) { rep = a.rep; return *this; }


// read-only access to representation
friend const GF2X& rep(const GF2E& a) { return a.rep; }

// You can always access the representation directly...if you dare.
GF2X& LoopHole() { return rep; }


static long WordLength() { return GF2EInfo->p.WordLength(); }

static const GF2XModulus& modulus() { return GF2EInfo->p; }

static long KarCross() { return GF2EInfo->KarCross; }
static long ModCross() { return GF2EInfo->ModCross; }
static long DivCross() { return GF2EInfo->DivCross; }

static long degree() { return GF2EInfo->p.n; }

static const GF2E& zero();

static const ZZ& cardinality() { return GF2EInfo->cardinality; }

static void init(const GF2X& NewP);

friend void clear(GF2E& x)
// x = 0
   { clear(x.rep); }

friend void set(GF2E& x)
// x = 1
   { set(x.rep); }

friend void swap(GF2E& x, GF2E& y)
// swap x and y

   { swap(x.rep, y.rep); }

// ****** addition

friend void add(GF2E& x, const GF2E& a, const GF2E& b)
   { add(x.rep, a.rep, b.rep); }

friend void add(GF2E& x, const GF2E& a, GF2 b)
   { add(x.rep, a.rep, b); }

friend void add(GF2E& x, const GF2E& a, long b)
   { add(x.rep, a.rep, b); }

friend void add(GF2E& x, GF2 a, const GF2E& b)  { add(x, b, a); }
friend void add(GF2E& x, long a, const GF2E& b)  { add(x, b, a); }

friend void sub(GF2E& x, const GF2E& a, const GF2E& b) { add(x, a, b); }
friend void sub(GF2E& x, const GF2E& a, GF2 b) { add(x, a, b); }
friend void sub(GF2E& x, const GF2E& a, long b) { add(x, a, b); }
friend void sub(GF2E& x, GF2 a, const GF2E& b) { add(x, a, b); }
friend void sub(GF2E& x, long a, const GF2E& b) { add(x, a, b); }

friend void negate(GF2E& x, const GF2E& a) { x = a; }


friend GF2E operator+(const GF2E& a, const GF2E& b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator+(const GF2E& a, GF2 b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator+(const GF2E& a, long b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator+(GF2 a, const GF2E& b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator+(long a, const GF2E& b) 
   { GF2E x; add(x, a, b); NTL_OPT_RETURN(GF2E, x); }


friend GF2E operator-(const GF2E& a, const GF2E& b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(const GF2E& a, GF2 b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(const GF2E& a, long b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(GF2 a, const GF2E& b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(long a, const GF2E& b) 
   { GF2E x; sub(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator-(const GF2E& a)
   { GF2E x; negate(x, a); NTL_OPT_RETURN(GF2E, x); } 


friend GF2E& operator+=(GF2E& x, const GF2E& b)
   { add(x, x, b); return x; }

friend GF2E& operator+=(GF2E& x, GF2 b)
   { add(x, x, b); return x; }

friend GF2E& operator+=(GF2E& x, long b)
   { add(x, x, b); return x; }


friend GF2E& operator-=(GF2E& x, const GF2E& b)
   { sub(x, x, b); return x; }

friend GF2E& operator-=(GF2E& x, GF2 b)
   { sub(x, x, b); return x; }

friend GF2E& operator-=(GF2E& x, long b)
   { sub(x, x, b); return x; }


friend GF2E& operator++(GF2E& x) { add(x, x, 1); return x; }

friend void operator++(GF2E& x, int) { add(x, x, 1); }

friend GF2E& operator--(GF2E& x) { sub(x, x, 1); return x; }

friend void operator--(GF2E& x, int) { sub(x, x, 1); }



// ****** multiplication

friend void mul(GF2E& x, const GF2E& a, const GF2E& b)
// x = a*b

   { MulMod(x.rep, a.rep, b.rep, GF2E::modulus()); }


friend void sqr(GF2E& x, const GF2E& a)
// x = a^2

   { SqrMod(x.rep, a.rep, GF2E::modulus()); }

friend GF2E sqr(const GF2E& a)
   { GF2E x; sqr(x, a); NTL_OPT_RETURN(GF2E, x); }

friend void mul(GF2E& x, const GF2E& a, GF2 b)
   { mul(x.rep, a.rep, b); }

friend void mul(GF2E& x, const GF2E& a, long b)
   { mul(x.rep, a.rep, b); }

friend void mul(GF2E& x, GF2 a, const GF2E& b) { mul(x, b, a); }
friend void mul(GF2E& x, long a, const GF2E& b) { mul(x, b, a); }



friend GF2E operator*(const GF2E& a, const GF2E& b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator*(const GF2E& a, GF2 b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator*(const GF2E& a, long b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator*(GF2 a, const GF2E& b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator*(long a, const GF2E& b) 
   { GF2E x; mul(x, a, b); NTL_OPT_RETURN(GF2E, x); }


friend GF2E& operator*=(GF2E& x, const GF2E& b)
   { mul(x, x, b); return x; }

friend GF2E& operator*=(GF2E& x, GF2 b)
   { mul(x, x, b); return x; }

friend GF2E& operator*=(GF2E& x, long b)
   { mul(x, x, b); return x; }



// ****** division



friend void div(GF2E& x, const GF2E& a, const GF2E& b);

friend void inv(GF2E& x, const GF2E& a);

friend GF2E inv(const GF2E& a)
   { GF2E x; inv(x, a); NTL_OPT_RETURN(GF2E, x); }

friend void div(GF2E& x, const GF2E& a, GF2 b)
   { div(x.rep, a.rep, b); } 

friend void div(GF2E& x, const GF2E& a, long b)
   { div(x.rep, a.rep, b); } 

friend void div(GF2E& x, GF2 a, const GF2E& b);
friend void div(GF2E& x, long a, const GF2E& b);


friend GF2E operator/(const GF2E& a, const GF2E& b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator/(const GF2E& a, GF2 b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator/(const GF2E& a, long b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator/(GF2 a, const GF2E& b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }

friend GF2E operator/(long a, const GF2E& b) 
   { GF2E x; div(x, a, b); NTL_OPT_RETURN(GF2E, x); }


friend GF2E& operator/=(GF2E& x, const GF2E& b)
   { div(x, x, b); return x; }

friend GF2E& operator/=(GF2E& x, GF2 b)
   { div(x, x, b); return x; }

friend GF2E& operator/=(GF2E& x, long b)
   { div(x, x, b); return x; }


// ****** exponentiation

friend void power(GF2E& x, const GF2E& a, const ZZ& e)
   { PowerMod(x.rep, a.rep, e, GF2E::modulus()); }

friend GF2E power(const GF2E& a, const ZZ& e)
   { GF2E x; power(x, a, e); NTL_OPT_RETURN(GF2E, x); }

friend void power(GF2E& x, const GF2E& a, long e)
   { PowerMod(x.rep, a.rep, e, GF2E::modulus()); }

friend GF2E power(const GF2E& a, long e)
   { GF2E x; power(x, a, e); NTL_OPT_RETURN(GF2E, x); }


// ****** conversion

friend void conv(GF2E& x, const GF2X& a)
// x = (a mod p)

   { rem(x.rep, a, GF2E::modulus()); }

friend void conv(GF2E& x, long a)
   { conv(x.rep, a); }

friend void conv(GF2E& x, GF2 a)
   { conv(x.rep, a); }

friend void conv(GF2E& x, const ZZ& a)
   { conv(x.rep, a); }

friend GF2E to_GF2E(const GF2X& a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

friend GF2E to_GF2E(long a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

friend GF2E to_GF2E(GF2 a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }

friend GF2E to_GF2E(const ZZ& a)
   { GF2E x; conv(x, a); NTL_OPT_RETURN(GF2E, x); }


// ****** comparison

friend long IsZero(const GF2E& a)
   { return IsZero(a.rep); }

friend long IsOne(const GF2E& a)
   { return IsOne(a.rep); }

friend long operator==(const GF2E& a, const GF2E& b)
   { return a.rep == b.rep; }

friend long operator==(const GF2E& a, GF2 b)
   { return a.rep == b; }

friend long operator==(const GF2E& a, long b)
   { return a.rep == b; }

friend long operator==(const GF2 a, const GF2E& b)
   { return a == b.rep; }

friend long operator==(const long a, const GF2E& b)
   { return a == b.rep; }


friend long operator!=(const GF2E& a, const GF2E& b) { return !(a == b); }
friend long operator!=(const GF2E& a, GF2 b) { return !(a == b); }
friend long operator!=(const GF2E& a, long b) { return !(a == b); }
friend long operator!=(GF2 a, const GF2E& b) { return !(a == b); }
friend long operator!=(long a, const GF2E& b) { return !(a == b); }

// ****** trace

friend void trace(GF2& x, const GF2E& a)
   { TraceMod(x, a.rep, GF2E::modulus()); }
friend GF2 trace(const GF2E& a)
   { return TraceMod(a.rep, GF2E::modulus()); }



// ****** random numbers

friend void random(GF2E& x)
// x = random element in GF2E

   { random(x.rep, GF2EInfo->p.n); }

friend GF2E random_GF2E()
   { GF2E x; random(x); NTL_OPT_RETURN(GF2E, x); }


// ****** input/output

friend ostream& operator<<(ostream& s, const GF2E& a)
   { return s << a.rep; }
   
friend istream& operator>>(istream& s, GF2E& x);


GF2E& operator=(long a) { conv(*this, a); return *this; }
GF2E& operator=(GF2 a) { conv(*this, a); return *this; }


// some other friends...not for general use

friend void BlockConstruct(GF2E*,long);
friend void BlockDestroy(GF2E*,long);


};


#endif
