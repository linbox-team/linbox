

#ifndef NTL_ZZ_p__H
#define NTL_ZZ_p__H

#include <NTL/ZZ.h>
#include <NTL/ZZVec.h>



// representation:  each ZZ_p is represented by a ZZ in the range 0..p-1.

// The constructor for a ZZ_p pre-allocates space for the underlying ZZ,
// and initializes it to zero.

const int MAX_ZZ_p_TEMPS = 16;

class ZZ_p;

class ZZ_pInfoT {
private:
   ZZ_pInfoT();                       // disabled
   ZZ_pInfoT(const ZZ_pInfoT&);   // disabled
   void operator=(const ZZ_pInfoT&);  // disabled
public:
   ZZ_pInfoT(const ZZ& NewP);
   ~ZZ_pInfoT();

   long ref_count; // reference count for gargabge collection

   ZZ p;      // the modulus
   long size;  // p.size()
   long ExtendedModulusSize;

   // the following implement a "lazy" initialization strategy
   long initialized;  // flag if initialization really was done
   void init();
   void check() { if (!initialized) init(); }

   long NumPrimes;

   long MaxRoot;  

   long QuickCRT;

   ZZ MinusMModP;  //  -M mod p, M = product of primes

   long mmsize;     // size pre-computed for MultiMul

   // the following arrays are indexed 0..NumPrimes-1
   // q = FFTPrime[i]


   ZZVec CoeffModP;    // coeff mod p

   double *x;          // u/q, where u = (M/q)^{-1} mod q
   long *u;            // u, as above

   double **tbl;       // table used for MultiRem; only with NTL_SINGLE_MUL

   long **tbl1;        // table used for MultiRem; only with NTL_TBL_REM

   ZZ_p *temps[MAX_ZZ_p_TEMPS];
   long temps_top;
};

extern ZZ_pInfoT *ZZ_pInfo; // info for current modulus, initially null



class ZZ_pContext {
private:
ZZ_pInfoT *ptr;

public:
void save();
void restore() const;

ZZ_pContext() { ptr = 0; }
ZZ_pContext(const ZZ& p);

ZZ_pContext(const ZZ_pContext&); 


ZZ_pContext& operator=(const ZZ_pContext&); 


~ZZ_pContext();


};


class ZZ_pBak {
private:
long MustRestore;
ZZ_pInfoT *ptr;

ZZ_pBak(const ZZ_pBak&); // disabled
void operator=(const ZZ_pBak&); // disabled

public:
void save();
void restore();

ZZ_pBak() { MustRestore = 0; ptr = 0; }

~ZZ_pBak();


};




struct ZZ_p_NoAlloc_type { ZZ_p_NoAlloc_type() { } };
const ZZ_p_NoAlloc_type ZZ_p_NoAlloc = ZZ_p_NoAlloc_type();


class ZZ_pTemp {
private:
   long pos;

public:
   ZZ_pTemp();
   ~ZZ_pTemp();

   ZZ_p& val() const;
};

#define NTL_ZZ_pRegister(x)  \
   ZZ_pTemp ZZ_pTemp__ ## x; ZZ_p& x = ZZ_pTemp__ ## x . val()


class ZZ_p {

private:

ZZ rep;

// static data

public:

static void init(const ZZ&);


typedef void (*DivHandlerPtr)(const ZZ_p& a);   // error-handler for division
static DivHandlerPtr DivHandler;


// ****** constructors and assignment

ZZ_p();

ZZ_p(const ZZ_p& a) :  rep(INIT_SIZE, ZZ_pInfo->size) { rep = a.rep; }

ZZ_p(ZZ_p_NoAlloc_type) { }  // allocates no space

~ZZ_p() { } 

ZZ_p& operator=(const ZZ_p& a) { rep = a.rep; return *this; }

// read-only access to representation
friend const ZZ& rep(const ZZ_p& a) { return a.rep; }

// You can always access the representation directly...if you dare.
ZZ& LoopHole() { return rep; }

ZZ_p(ZZ_p& x, INIT_TRANS_TYPE) : rep(x.rep, INIT_TRANS) { }



static const ZZ& modulus() { return ZZ_pInfo->p; }
static long ModulusSize() { return ZZ_pInfo->size; }
static long StorageSize() { return ZZ_pInfo->p.size()+4; }

static const ZZ_p& zero();

// ****** conversion

friend void conv(ZZ_p& x, const ZZ& a)
   { rem(x.rep, a, ZZ_p::modulus()); }

ZZ_p(INIT_VAL_TYPE, const ZZ& a);

friend ZZ_p to_ZZ_p(const ZZ& a)
   { return ZZ_p(INIT_VAL, a); }


friend void conv(ZZ_p& x, long a);

ZZ_p(INIT_VAL_TYPE, long a);

friend ZZ_p to_ZZ_p(long a)
   { return ZZ_p(INIT_VAL, a); }


// ****** assignment

ZZ_p& operator=(long a) { conv(*this, a); return *this; }


// ****** some basics


friend void clear(ZZ_p& x)
// x = 0
   { clear(x.rep); }

friend void set(ZZ_p& x)
// x = 1
   { set(x.rep); }

friend void swap(ZZ_p& x, ZZ_p& y)
// swap x and y

   { swap(x.rep, y.rep); }

// ****** addition

friend void add(ZZ_p& x, const ZZ_p& a, const ZZ_p& b)
// x = a + b

   { AddMod(x.rep, a.rep, b.rep, ZZ_p::modulus()); }

friend void sub(ZZ_p& x, const ZZ_p& a, const ZZ_p& b)
// x = a - b

   { SubMod(x.rep, a.rep, b.rep, ZZ_p::modulus()); }


friend void negate(ZZ_p& x, const ZZ_p& a)
// x = -a

   { NegateMod(x.rep, a.rep, ZZ_p::modulus()); }


// scalar versions

friend void add(ZZ_p& x, const ZZ_p& a, long b);
friend void add(ZZ_p& x, long a, const ZZ_p& b) { add(x, b, a); }

friend void sub(ZZ_p& x, const ZZ_p& a, long b);
friend void sub(ZZ_p& x, long a, const ZZ_p& b);


// ****** multiplication

friend void mul(ZZ_p& x, const ZZ_p& a, const ZZ_p& b)
// x = a*b

   { MulMod(x.rep, a.rep, b.rep, ZZ_p::modulus()); }


friend void sqr(ZZ_p& x, const ZZ_p& a)
// x = a^2

   { SqrMod(x.rep, a.rep, ZZ_p::modulus()); }

friend ZZ_p sqr(const ZZ_p& a)
   { ZZ_p x; sqr(x, a); NTL_OPT_RETURN(ZZ_p, x); }


// scalar versions

friend void mul(ZZ_p& x, const ZZ_p& a, long b);
friend void mul(ZZ_p& x, long a, const ZZ_p& b) { mul(x, b, a); }

// ****** division


friend void div(ZZ_p& x, const ZZ_p& a, const ZZ_p& b);
// x = a/b
// If b != 0 & b not invertible & DivHandler != 0,
// then DivHandler will be called with the offending b.
// In this case, of course, p is not really prime, and one
// can factor p by taking a gcd with rep(b).
// Otherwise, if b is not invertible, an error occurs.

friend void inv(ZZ_p& x, const ZZ_p& a);
// x = 1/a
// Error handling is the same as above.

friend ZZ_p inv(const ZZ_p& a)
   { ZZ_p x; inv(x, a); NTL_OPT_RETURN(ZZ_p, x); }

friend void div(ZZ_p& x, const ZZ_p& a, long b);
friend void div(ZZ_p& x, long a, const ZZ_p& b);


// operator notation:

friend ZZ_p operator+(const ZZ_p& a, const ZZ_p& b)
   { ZZ_p x; add(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p operator+(const ZZ_p& a, long b)
   { ZZ_p x; add(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p operator+(long a, const ZZ_p& b)
   { ZZ_p x; add(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p& operator+=(ZZ_p& x, const ZZ_p& b)
   { add(x, x, b); return x; } 

friend ZZ_p& operator+=(ZZ_p& x, long b)
   { add(x, x, b); return x; } 



friend ZZ_p operator-(const ZZ_p& a, const ZZ_p& b)
   { ZZ_p x; sub(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p operator-(const ZZ_p& a, long b)
   { ZZ_p x; sub(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p operator-(long a, const ZZ_p& b)
   { ZZ_p x; sub(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p& operator-=(ZZ_p& x, const ZZ_p& b)
   { sub(x, x, b); return x; } 

friend ZZ_p& operator-=(ZZ_p& x, long b)
   { sub(x, x, b); return x; } 



friend ZZ_p operator*(const ZZ_p& a, const ZZ_p& b)
   { ZZ_p x; mul(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p operator*(const ZZ_p& a, long b)
   { ZZ_p x; mul(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p operator*(long a, const ZZ_p& b)
   { ZZ_p x; mul(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p& operator*=(ZZ_p& x, const ZZ_p& b)
   { mul(x, x, b); return x; } 

friend ZZ_p& operator*=(ZZ_p& x, long b)
   { mul(x, x, b); return x; } 


friend ZZ_p operator/(const ZZ_p& a, const ZZ_p& b)
   { ZZ_p x; div(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p operator/(const ZZ_p& a, long b)
   { ZZ_p x; div(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p operator/(long a, const ZZ_p& b)
   { ZZ_p x; div(x, a, b); NTL_OPT_RETURN(ZZ_p, x); }

friend ZZ_p& operator/=(ZZ_p& x, const ZZ_p& b)
   { div(x, x, b); return x; } 

friend ZZ_p& operator/=(ZZ_p& x, long b)
   { div(x, x, b); return x; } 


friend ZZ_p operator-(const ZZ_p& a)
   { ZZ_p x; negate(x, a); NTL_OPT_RETURN(ZZ_p, x); }


friend ZZ_p& operator++(ZZ_p& x) { add(x, x, 1); return x; }
friend void operator++(ZZ_p& x, int) { add(x, x, 1); }
friend ZZ_p& operator--(ZZ_p& x) { sub(x, x, 1); return x; }
friend void operator--(ZZ_p& x, int) { sub(x, x, 1); }


// ****** exponentiation

friend void power(ZZ_p& x, const ZZ_p& a, const ZZ& e)
   { PowerMod(x.rep, a.rep, e, ZZ_p::modulus()); }

friend ZZ_p power(const ZZ_p& a, const ZZ& e)
   { ZZ_p x; power(x, a, e); NTL_OPT_RETURN(ZZ_p, x); }

friend void power(ZZ_p& x, const ZZ_p& a, long e)
   { PowerMod(x.rep, a.rep, e, ZZ_p::modulus()); }

friend ZZ_p power(const ZZ_p& a, long e)
   { ZZ_p x; power(x, a, e); NTL_OPT_RETURN(ZZ_p, x); }


// ****** comparison

friend long IsZero(const ZZ_p& a)
   { return IsZero(a.rep); }


friend long IsOne(const ZZ_p& a)
   { return IsOne(a.rep); }

friend long operator==(const ZZ_p& a, const ZZ_p& b)
   { return a.rep == b.rep; }

friend long operator!=(const ZZ_p& a, const ZZ_p& b)
   { return !(a == b); }

friend long operator==(const ZZ_p& a, long b);
friend long operator==(long a, const ZZ_p& b) { return b == a; }

friend long operator!=(const ZZ_p& a, long b) { return !(a == b); }
friend long operator!=(long a, const ZZ_p& b) { return !(a == b); }


// ****** random numbers

friend void random(ZZ_p& x)
// x = random element in ZZ_p

   { RandomBnd(x.rep, ZZ_p::modulus()); }

friend ZZ_p random_ZZ_p()
   { ZZ_p x; random(x); NTL_OPT_RETURN(ZZ_p, x); }


// ****** input/output

friend ostream& operator<<(ostream& s, const ZZ_p& a)
   { return s << a.rep; }
   
friend istream& operator>>(istream& s, ZZ_p& x);


// some other friends...not for general use

friend void BlockConstruct(ZZ_p*,long);
friend void BlockDestroy(ZZ_p*,long);


};


#endif
