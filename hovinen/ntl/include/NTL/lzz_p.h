
#ifndef NTL_zz_p__H
#define NTL_zz_p__H

#include <NTL/ZZ.h>
#include <NTL/FFT.h>



class zz_pInfoT {
private:
   zz_pInfoT();                      // disabled
   zz_pInfoT(const zz_pInfoT&);  // disabled
   void operator=(const zz_pInfoT&); // disabled
public:
   zz_pInfoT(long NewP, long maxroot);
   zz_pInfoT(long Index);
   ~zz_pInfoT();

   long ref_count;

   long p;
   double pinv;

   long index;        // index >= 0 means we are directly using
                     // an FFT prime

   long PrimeCnt;     // 0 for FFT prime;  otherwise same as NumPrimes
                     // used for establishing crossover points

   long NumPrimes;

   long MaxRoot;

   long MinusMModP;  //  -M mod p, M = product of primes

   // the following arrays are indexed 0..NumPrimes-1
   // q = FFTPrime[i]


   long *CoeffModP;    // coeff mod p

   double *x;          // u/q, where u = (M/q)^{-1} mod q
   long *u;            // u, as above
};

extern zz_pInfoT *zz_pInfo;  // current modulus, initially null



class zz_pContext {
private:
zz_pInfoT *ptr;

public:
void save();
void restore() const;

zz_pContext() { ptr = 0; }
zz_pContext(long p, long maxroot=NTL_FFTMaxRoot);
zz_pContext(INIT_FFT_TYPE, long index);

zz_pContext(const zz_pContext&); 

zz_pContext& operator=(const zz_pContext&); 

~zz_pContext();


};


class zz_pBak {
private:
long MustRestore;
zz_pInfoT *ptr;

zz_pBak(const zz_pBak&); // disabled
void operator=(const zz_pBak&); // disabled

public:
void save();
void restore();

zz_pBak() { MustRestore = 0; ptr = 0; }

~zz_pBak();


};

#define NTL_zz_pRegister(x) zz_p x


class zz_p {

long rep;

// static data


public:

static void init(long NewP, long maxroot=NTL_FFTMaxRoot);
static void FFTInit(long index);

friend class zz_pBak;


// ****** constructors and assignment

zz_p() { rep = 0; }

zz_p(const zz_p& a) :  rep(a.rep) { }  

~zz_p() { } 

zz_p& operator=(const zz_p& a) { rep = a.rep; return *this; }

friend zz_p to_zz_p(long a);
friend void conv(zz_p& x, long a);

friend zz_p to_zz_p(const ZZ& a);
friend void conv(zz_p& x, const ZZ& a);

zz_p& operator=(long a) { conv(*this, a); return *this; }

// read-only access to representation
friend long rep(zz_p a) { return a.rep; }

// a loop-hole for direct access to rep
long& LoopHole() { return rep; }

static long modulus() { return zz_pInfo->p; }
static zz_p zero() { return zz_p(); }
static double ModulusInverse() { return zz_pInfo->pinv; }
static long PrimeCnt() { return zz_pInfo->PrimeCnt; }

static long StorageSize() { return 1; }

friend void clear(zz_p& x)
// x = 0
   { x.rep = 0; }

friend void set(zz_p& x)
// x = 1
   { x.rep = 1; }

friend void swap(zz_p& x, zz_p& y)
// swap x and y

   { long t;  t = x.rep; x.rep = y.rep; y.rep = t; }

// ****** addition

friend void add(zz_p& x, zz_p a, zz_p b)
// x = a + b

   { x.rep = AddMod(a.rep, b.rep, zz_p::modulus()); }

friend void sub(zz_p& x, zz_p a, zz_p b)
// x = a - b

   { x.rep = SubMod(a.rep, b.rep, zz_p::modulus()); }


friend void negate(zz_p& x, zz_p a)
// x = -a

   { x.rep = SubMod(0, a.rep, zz_p::modulus()); }

// scalar versions

friend void add(zz_p& x, zz_p a, long b) { add(x, a, to_zz_p(b)); }
friend void add(zz_p& x, long a, zz_p b) { add(x, to_zz_p(a), b); }

friend void sub(zz_p& x, zz_p a, long b) { sub(x, a, to_zz_p(b)); }
friend void sub(zz_p& x, long a, zz_p b) { sub(x, to_zz_p(a), b); }

friend zz_p operator+(zz_p a, zz_p b)
    { zz_p x; add(x, a, b); return x; }

friend zz_p operator+(zz_p a, long b)
    { zz_p x; add(x, a, b); return x; }

friend zz_p operator+(long a, zz_p b)
    { zz_p x; add(x, a, b); return x; }

friend zz_p operator-(zz_p a, zz_p b)
    { zz_p x; sub(x, a, b); return x; }

friend zz_p operator-(zz_p a, long b)
    { zz_p x; sub(x, a, b); return x; }

friend zz_p operator-(long a, zz_p b)
    { zz_p x; sub(x, a, b); return x; }



friend zz_p operator-(zz_p a)
   { zz_p x; negate(x, a); return x; }



friend zz_p& operator+=(zz_p& x, zz_p b)
   { add(x, x, b); return x; }

friend zz_p& operator+=(zz_p& x, long b)
   { add(x, x, b); return x; }



friend zz_p& operator-=(zz_p& x, zz_p b)
   { sub(x, x, b); return x; }

friend zz_p& operator-=(zz_p& x, long b)
   { sub(x, x, b); return x; }

friend zz_p& operator++(zz_p& x) { add(x, x, 1); return x; }
friend void operator++(zz_p& x, int) { add(x, x, 1); }
friend zz_p& operator--(zz_p& x) { sub(x, x, 1); return x; }
friend void operator--(zz_p& x, int) { sub(x, x, 1); }

// ****** multiplication

friend void mul(zz_p& x, zz_p a, zz_p b)
// x = a*b

   { x.rep = MulMod(a.rep, b.rep, zz_p::modulus(), zz_p::ModulusInverse()); }

friend void mul(zz_p& x, zz_p a, long b) { mul(x, a, to_zz_p(b)); }
friend void mul(zz_p& x, long a, zz_p b) { mul(x, to_zz_p(a), b); }

friend zz_p operator*(zz_p a, zz_p b)
    { zz_p x; mul(x, a, b); return x; }

friend zz_p operator*(zz_p a, long b)
    { zz_p x; mul(x, a, b); return x; }

friend zz_p operator*(long a, zz_p b)
    { zz_p x; mul(x, a, b); return x; }


friend zz_p& operator*=(zz_p& x, zz_p b)
   { mul(x, x, b); return x; }

friend zz_p& operator*=(zz_p& x, long b)
   { mul(x, x, b); return x; }



friend void sqr(zz_p& x, zz_p a)
// x = a^2

   { x.rep = MulMod(a.rep, a.rep, zz_p::modulus(), zz_p::ModulusInverse()); }

friend zz_p sqr(zz_p a)
   { zz_p x; sqr(x, a); return x; }



// ****** division

friend void div(zz_p& x, zz_p a, zz_p b)
// x = a/b

   { x.rep = MulMod(a.rep, InvMod(b.rep, zz_p::modulus()), zz_p::modulus(),
                    zz_p::ModulusInverse()); }

friend void inv(zz_p& x, zz_p a)
// x = 1/a

   { x.rep = InvMod(a.rep, zz_p::modulus()); }

friend zz_p inv(zz_p a)
   { zz_p x; inv(x, a); return x; }

friend void div(zz_p& x, zz_p a, long b) { div(x, a, to_zz_p(b)); }
friend void div(zz_p& x, long a, zz_p b) { div(x, to_zz_p(a), b); }

friend zz_p operator/(zz_p a, zz_p b)
    { zz_p x; div(x, a, b); return x; }

friend zz_p operator/(zz_p a, long b)
    { zz_p x; div(x, a, b); return x; }

friend zz_p operator/(long a, zz_p b)
    { zz_p x; div(x, a, b); return x; }


friend zz_p& operator/=(zz_p& x, zz_p b)
   { div(x, x, b); return x; }

friend zz_p& operator/=(zz_p& x, long b)
   { div(x, x, b); return x; }


// ****** exponentiation

friend void power(zz_p& x, zz_p a, long e)
// x = a^e

   { x.rep = PowerMod(a.rep, e, zz_p::modulus()); }

friend zz_p power(zz_p a, long e)
   { zz_p x; power(x, a, e); return x; }

// ****** comparison

friend long IsZero(zz_p a)
   { return a.rep == 0; }

friend long IsOne(zz_p a)
   { return a.rep == 1; }

friend long operator==(zz_p a, zz_p b)
   { return a.rep == b.rep; }

friend long operator!=(zz_p a, zz_p b)
   { return !(a == b); }

friend long operator==(zz_p a, long b) { return a == to_zz_p(b); }
friend long operator==(long a, zz_p b) { return to_zz_p(a) == b; }

friend long operator!=(zz_p a, long b) { return !(a == b); }
friend long operator!=(long a, zz_p b) { return !(a == b); }

// ****** random numbers

friend void random(zz_p& x)
// x = random element in zz_p

   { x.rep = RandomBnd(zz_p::modulus()); }

friend zz_p random_zz_p()
   { zz_p x; random(x); return x; }


// ****** input/output

friend ostream& operator<<(ostream& s, zz_p a);
   
friend istream& operator>>(istream& s, zz_p& x);

zz_p(long a, INIT_LOOP_HOLE_TYPE) { rep = a; }

};


#endif
