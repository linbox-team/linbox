
#ifndef NTL_ZZ_pE__H
#define NTL_ZZ_pE__H

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/vec_long.h>

#include <NTL/ZZ_pX.h>


class ZZ_pE;

class ZZ_pEInfoT {
private:
   ZZ_pEInfoT();                       // disabled
   ZZ_pEInfoT(const ZZ_pEInfoT&);   // disabled
   void operator=(const ZZ_pEInfoT&);  // disabled
public:
   long ref_count;

   ZZ_pEInfoT(const ZZ_pX&);
   ~ZZ_pEInfoT() { }

   ZZ_pXModulus p;

   ZZ cardinality;

void add(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b);
void sub(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b);
void negate(ZZ_pE& x, const ZZ_pE& a);

void add(ZZ_pE& x, const ZZ_pE& a, long b);
void add(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b);
void add(ZZ_pE& x, long a, const ZZ_pE& b);
void add(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b);

void sub(ZZ_pE& x, const ZZ_pE& a, long b);
void sub(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b);
void sub(ZZ_pE& x, long a, const ZZ_pE& b);
void sub(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b);

void mul(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b);
void sqr(ZZ_pE& x, const ZZ_pE& a);
ZZ_pE sqr(const ZZ_pE& a);

void mul(ZZ_pE& x, const ZZ_pE& a, long b);
void mul(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b);
void mul(ZZ_pE& x, long a, const ZZ_pE& b);
void mul(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b);

void div(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b);
void div(ZZ_pE& x, const ZZ_pE& a, long b);
void div(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b);
void div(ZZ_pE& x, long a, const ZZ_pE& b);
void div(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b);

void inv(ZZ_pE& x, const ZZ_pE& a);
ZZ_pE inv(const ZZ_pE& a);

void power(ZZ_pE& x, const ZZ_pE& a, const ZZ& e);
ZZ_pE power(const ZZ_pE& a, const ZZ& e);
void power(ZZ_pE& x, const ZZ_pE& a, long e);
ZZ_pE power(const ZZ_pE& a, long e);

void conv(ZZ_pE& x, const ZZ_pX& a);
void conv(ZZ_pE& x, long a);
void conv(ZZ_pE& x, const ZZ_p& a);
void conv(ZZ_pE& x, const ZZ& a);
ZZ_pE to_ZZ_pE(const ZZ_pX& a);
ZZ_pE to_ZZ_pE(long a);
ZZ_pE to_ZZ_pE(const ZZ_p& a);
ZZ_pE to_ZZ_pE(const ZZ& a);

void trace(ZZ_p& x, const ZZ_pE& a);
ZZ_p trace(const ZZ_pE& a);

void norm(ZZ_p& x, const ZZ_pE& a);
ZZ_p norm(const ZZ_pE& a);

void random(ZZ_pE& x);
ZZ_pE random_ZZ_pE();

};

extern ZZ_pEInfoT *ZZ_pEInfo; // info for current modulus, initially null


#if (defined (_THREAD_SAFE)) || (defined (_REENTRANT))
#  if !defined (COARSE_LOCKS)

extern pthread_rwlock_t ZZ_pE_lock;

#  else
#    define ZZ_pE_lock field_lock;
#  endif
#  define NTL_ZZ_PE_THREADS_ENTER pthread_rwlock_rdlock (&ZZ_pE_lock);
#  define NTL_ZZ_PE_THREADS_LEAVE pthread_rwlock_unlock (&ZZ_pE_lock);
#else
#  define NTL_ZZ_PE_THREADS_ENTER
#  define NTL_ZZ_PE_THREADS_LEAVE
#endif



class ZZ_pEContext {
private:
ZZ_pEInfoT *ptr;

public:
void save();
void restore() const;

ZZ_pEContext() { ptr = 0; }
ZZ_pEContext(const ZZ_pX& p);

ZZ_pEContext(const ZZ_pEContext&); 


ZZ_pEContext& operator=(const ZZ_pEContext&); 


~ZZ_pEContext();

};


class ZZ_pEBak {
private:
long MustRestore;
ZZ_pEInfoT *ptr;

ZZ_pEBak(const ZZ_pEBak&); // disabled
void operator=(const ZZ_pEBak&); // disabled

public:
void save();
void restore();

ZZ_pEBak() { MustRestore = 0; ptr = 0; }

~ZZ_pEBak();


};



struct ZZ_pE_NoAlloc_type { ZZ_pE_NoAlloc_type() { } };
const ZZ_pE_NoAlloc_type ZZ_pE_NoAlloc = ZZ_pE_NoAlloc_type();



class ZZ_pE {


private:

friend class ZZ_pEInfoT;

ZZ_pX rep;

// static data


public:

static long DivCross() { return 16; }
static long ModCross() { return 8; }


// ****** constructors and assignment

ZZ_pE();

ZZ_pE(const ZZ_pE& a)  { rep.rep.SetMaxLength(ZZ_pE::degree()); rep = a.rep; }

ZZ_pE(ZZ_pE_NoAlloc_type) { }  // allocates no space

~ZZ_pE() { } 

ZZ_pE& operator=(const ZZ_pE& a) { rep = a.rep; return *this; }

ZZ_pE(ZZ_pE& x, INIT_TRANS_TYPE) : rep(x.rep, INIT_TRANS) { }


// read-only access to representation
friend const ZZ_pX& rep(const ZZ_pE& a) { return a.rep; }

// You can always access the representation directly...if you dare.
ZZ_pX& LoopHole() { return rep; }



static const ZZ_pXModulus& modulus() { return ZZ_pEInfo->p; }

static long degree() { return deg(ZZ_pEInfo->p); }

static const ZZ& cardinality() { return ZZ_pEInfo->cardinality; }

static const ZZ_pE& zero();

static long initialized() { return (ZZ_pEInfo != 0); }

static void init(const ZZ_pX&);

friend void clear(ZZ_pE& x)
// x = 0
   { clear(x.rep); }

friend void set(ZZ_pE& x)
// x = 1
   { set(x.rep); }

friend void swap(ZZ_pE& x, ZZ_pE& y)
// swap x and y

   { swap(x.rep, y.rep); }

// ****** addition

friend void add(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a + b

   { add(x.rep, a.rep, b.rep); }

friend void sub(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a - b

   { sub(x.rep, a.rep, b.rep); }


friend void negate(ZZ_pE& x, const ZZ_pE& a) 

   { negate(x.rep, a.rep); }


friend void add(ZZ_pE& x, const ZZ_pE& a, long b)
   { add(x.rep, a.rep, b); }

friend void add(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { add(x.rep, a.rep, b); }

friend void add(ZZ_pE& x, long a, const ZZ_pE& b)
   { add(x.rep, a, b.rep); }

friend void add(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { add(x.rep, a, b.rep); }





friend void sub(ZZ_pE& x, const ZZ_pE& a, long b)
   { sub(x.rep, a.rep, b); }

friend void sub(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { sub(x.rep, a.rep, b); }

friend void sub(ZZ_pE& x, long a, const ZZ_pE& b)
   { sub(x.rep, a, b.rep); }

friend void sub(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { sub(x.rep, a, b.rep); }





// ****** multiplication

friend void mul(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a*b

   { MulMod(x.rep, a.rep, b.rep, ZZ_pE::modulus()); }


friend void sqr(ZZ_pE& x, const ZZ_pE& a)
// x = a^2

   { SqrMod(x.rep, a.rep, ZZ_pE::modulus()); }

friend ZZ_pE sqr(const ZZ_pE& a)
   { ZZ_pE x; sqr(x, a); NTL_OPT_RETURN(ZZ_pE, x); }


friend void mul(ZZ_pE& x, const ZZ_pE& a, long b)
   { mul(x.rep, a.rep, b); }

friend void mul(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { mul(x.rep, a.rep, b); }

friend void mul(ZZ_pE& x, long a, const ZZ_pE& b)
   { mul(x.rep, a, b.rep); }

friend void mul(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { mul(x.rep, a, b.rep); }


// ****** division




friend void inv(ZZ_pE& x, const ZZ_pE& a);
friend ZZ_pE inv(const ZZ_pE& a)
   { ZZ_pE x; inv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }



// ****** exponentiation

friend void power(ZZ_pE& x, const ZZ_pE& a, const ZZ& e)
// x = a^e

   { PowerMod(x.rep, a.rep, e, ZZ_pE::modulus()); }

friend ZZ_pE power(const ZZ_pE& a, const ZZ& e)
   { ZZ_pE x; power(x, a, e); NTL_OPT_RETURN(ZZ_pE, x); }

friend void power(ZZ_pE& x, const ZZ_pE& a, long e)
   { _BUFFER ZZ ee; conv (ee, e); power(x, a, ee); }

friend ZZ_pE power(const ZZ_pE& a, long e)
   { ZZ_pE x; power(x, a, e); NTL_OPT_RETURN(ZZ_pE, x); }




// ****** conversion

friend void conv(ZZ_pE& x, const ZZ_pX& a)
   { rem(x.rep, a, ZZ_pE::modulus()); }

friend void conv(ZZ_pE& x, long a)
   { conv(x.rep, a); }

friend void conv(ZZ_pE& x, const ZZ_p& a)
   { conv(x.rep, a); }

friend void conv(ZZ_pE& x, const ZZ& a)
   { conv(x.rep, a); }

friend ZZ_pE to_ZZ_pE(const ZZ_pX& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

friend ZZ_pE to_ZZ_pE(long a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

friend ZZ_pE to_ZZ_pE(const ZZ_p& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

friend ZZ_pE to_ZZ_pE(const ZZ& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }



// ****** comparison

friend long IsZero(const ZZ_pE& a)
   { return IsZero(a.rep); }

friend long IsOne(const ZZ_pE& a)
   { return IsOne(a.rep); }

friend long operator==(const ZZ_pE& a, const ZZ_pE& b)
   { return a.rep == b.rep; }
friend long operator==(const ZZ_pE& a, long b)
   { return a.rep == b; }
friend long operator==(const ZZ_pE& a, const ZZ_p& b)
   { return a.rep == b; }
friend long operator==(long a, const ZZ_pE& b)
   { return a == b.rep; }
friend long operator==(const ZZ_p& a, const ZZ_pE& b)
   { return a == b.rep; }

friend long operator!=(const ZZ_pE& a, const ZZ_pE& b)
   { return !(a == b); }
friend long operator!=(const ZZ_pE& a, long b)
   { return !(a == b); }
friend long operator!=(const ZZ_pE& a, const ZZ_p& b)
   { return !(a == b); }
friend long operator!=(long a, const ZZ_pE& b)
   { return !(a == b); }
friend long operator!=(const ZZ_p& a, const ZZ_pE& b)
   { return !(a == b); }


// ****** norm and trace

friend void trace(ZZ_p& x, const ZZ_pE& a)
   { TraceMod(x, a.rep, ZZ_pE::modulus()); }
friend ZZ_p trace(const ZZ_pE& a)
   { return TraceMod(a.rep, ZZ_pE::modulus()); }

friend void norm(ZZ_p& x, const ZZ_pE& a)
   { NormMod(x, a.rep, ZZ_pE::modulus()); }
friend ZZ_p norm(const ZZ_pE& a)
   { return NormMod(a.rep, ZZ_pE::modulus()); }


// ****** random numbers

friend void random(ZZ_pE& x)
// x = random element in ZZ_pE

   { random(x.rep, ZZ_pE::degree()); }

friend ZZ_pE random_ZZ_pE()
   { ZZ_pE x; random(x); NTL_OPT_RETURN(ZZ_pE, x); }


// ****** input/output

friend ostream& operator<<(ostream& s, const ZZ_pE& a)
   { return s << a.rep; }
   
friend istream& operator>>(istream& s, ZZ_pE& x);


ZZ_pE& operator=(long a) { conv(*this, a); return *this; }
ZZ_pE& operator=(const ZZ_p& a) { conv(*this, a); return *this; }


};

inline void div(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b) { ZZ_pEInfo->div (x, a, b); }
inline void div(ZZ_pE& x, const ZZ_pE& a, long b) { ZZ_pEInfo->div (x, a, b); }
inline void div(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b) { ZZ_pEInfo->div (x, a, b); }
inline void div(ZZ_pE& x, long a, const ZZ_pE& b) { ZZ_pEInfo->div (x, a, b); }
inline void div(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b) { ZZ_pEInfo->div (x, a, b); }

inline void inv(ZZ_pE& x, const ZZ_pE& a) { ZZ_pEInfo->inv (x, a); }


inline ZZ_pE operator+(const ZZ_pE& a, const ZZ_pE& b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator+(const ZZ_pE& a, const ZZ_p& b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator+(const ZZ_pE& a, long b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator+(const ZZ_p& a, const ZZ_pE& b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator+(long a, const ZZ_pE& b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }


inline ZZ_pE operator-(const ZZ_pE& a, const ZZ_pE& b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(const ZZ_pE& a, const ZZ_p& b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(const ZZ_pE& a, long b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(const ZZ_p& a, const ZZ_pE& b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(long a, const ZZ_pE& b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(const ZZ_pE& a)
   { ZZ_pE x; negate(x, a); NTL_OPT_RETURN(ZZ_pE, x); } 


inline ZZ_pE& operator+=(ZZ_pE& x, const ZZ_pE& b)
   { add(x, x, b); return x; }

inline ZZ_pE& operator+=(ZZ_pE& x, const ZZ_p& b)
   { add(x, x, b); return x; }

inline ZZ_pE& operator+=(ZZ_pE& x, long b)
   { add(x, x, b); return x; }


inline ZZ_pE& operator-=(ZZ_pE& x, const ZZ_pE& b)
   { sub(x, x, b); return x; }

inline ZZ_pE& operator-=(ZZ_pE& x, const ZZ_p& b)
   { sub(x, x, b); return x; }

inline ZZ_pE& operator-=(ZZ_pE& x, long b)
   { sub(x, x, b); return x; }


inline ZZ_pE& operator++(ZZ_pE& x) { add(x, x, 1); return x; }

inline void operator++(ZZ_pE& x, int) { add(x, x, 1); }

inline ZZ_pE& operator--(ZZ_pE& x) { sub(x, x, 1); return x; }

inline void operator--(ZZ_pE& x, int) { sub(x, x, 1); }



inline ZZ_pE operator*(const ZZ_pE& a, const ZZ_pE& b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator*(const ZZ_pE& a, const ZZ_p& b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator*(const ZZ_pE& a, long b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator*(const ZZ_p& a, const ZZ_pE& b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator*(long a, const ZZ_pE& b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }


inline ZZ_pE& operator*=(ZZ_pE& x, const ZZ_pE& b)
   { mul(x, x, b); return x; }

inline ZZ_pE& operator*=(ZZ_pE& x, const ZZ_p& b)
   { mul(x, x, b); return x; }

inline ZZ_pE& operator*=(ZZ_pE& x, long b)
   { mul(x, x, b); return x; }




inline ZZ_pE operator/(const ZZ_pE& a, const ZZ_pE& b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator/(const ZZ_pE& a, const ZZ_p& b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator/(const ZZ_pE& a, long b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator/(const ZZ_p& a, const ZZ_pE& b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator/(long a, const ZZ_pE& b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }


inline ZZ_pE& operator/=(ZZ_pE& x, const ZZ_pE& b)
   { div(x, x, b); return x; }

inline ZZ_pE& operator/=(ZZ_pE& x, const ZZ_p& b)
   { div(x, x, b); return x; }

inline ZZ_pE& operator/=(ZZ_pE& x, long b)
   { div(x, x, b); return x; }





// Versions of the above operations that are methods of ZZ_pEInfoT




inline void ZZ_pEInfoT::add(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a + b

   { ::add(x.rep, a.rep, b.rep); }

inline void ZZ_pEInfoT::sub(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a - b

   { ::sub(x.rep, a.rep, b.rep); }


inline void ZZ_pEInfoT::negate(ZZ_pE& x, const ZZ_pE& a) 

   { ::negate(x.rep, a.rep); }


inline void ZZ_pEInfoT::add(ZZ_pE& x, const ZZ_pE& a, long b)
   { ::add(x.rep, a.rep, b); }

inline void ZZ_pEInfoT::add(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { ::add(x.rep, a.rep, b); }

inline void ZZ_pEInfoT::add(ZZ_pE& x, long a, const ZZ_pE& b)
   { ::add(x.rep, a, b.rep); }

inline void ZZ_pEInfoT::add(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { ::add(x.rep, a, b.rep); }





inline void ZZ_pEInfoT::sub(ZZ_pE& x, const ZZ_pE& a, long b)
   { ::sub(x.rep, a.rep, b); }

inline void ZZ_pEInfoT::sub(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { ::sub(x.rep, a.rep, b); }

inline void ZZ_pEInfoT::sub(ZZ_pE& x, long a, const ZZ_pE& b)
   { ::sub(x.rep, a, b.rep); }

inline void ZZ_pEInfoT::sub(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { ::sub(x.rep, a, b.rep); }





// ****** multiplication

inline void ZZ_pEInfoT::mul(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a*b

   { MulMod(x.rep, a.rep, b.rep, p); }


inline void ZZ_pEInfoT::sqr(ZZ_pE& x, const ZZ_pE& a)
// x = a^2

   { SqrMod(x.rep, a.rep, p); }

inline ZZ_pE ZZ_pEInfoT::sqr(const ZZ_pE& a)
   { ZZ_pE x; sqr(x, a); NTL_OPT_RETURN(ZZ_pE, x); }


inline void ZZ_pEInfoT::mul(ZZ_pE& x, const ZZ_pE& a, long b)
   { ::mul(x.rep, a.rep, b); }

inline void ZZ_pEInfoT::mul(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { ::mul(x.rep, a.rep, b); }

inline void ZZ_pEInfoT::mul(ZZ_pE& x, long a, const ZZ_pE& b)
   { ::mul(x.rep, a, b.rep); }

inline void ZZ_pEInfoT::mul(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { ::mul(x.rep, a, b.rep); }

// ****** division

inline ZZ_pE ZZ_pEInfoT::inv(const ZZ_pE& a)
   { ZZ_pE x; inv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

// ****** exponentiation

inline void ZZ_pEInfoT::power(ZZ_pE& x, const ZZ_pE& a, const ZZ& e)
// x = a^e

   { PowerMod(x.rep, a.rep, e, p); }

inline ZZ_pE ZZ_pEInfoT::power(const ZZ_pE& a, const ZZ& e)
   { ZZ_pE x; power(x, a, e); NTL_OPT_RETURN(ZZ_pE, x); }

inline void ZZ_pEInfoT::power(ZZ_pE& x, const ZZ_pE& a, long e)
   { _BUFFER ZZ ee; ::conv (ee, e); power(x, a, ee); }

inline ZZ_pE ZZ_pEInfoT::power(const ZZ_pE& a, long e)
   { ZZ_pE x; power(x, a, e); NTL_OPT_RETURN(ZZ_pE, x); }




// ****** conversion

inline void ZZ_pEInfoT::conv(ZZ_pE& x, const ZZ_pX& a)
   { ::rem(x.rep, a, p); }

inline void ZZ_pEInfoT::conv(ZZ_pE& x, long a)
   { ::conv(x.rep, a); }

inline void ZZ_pEInfoT::conv(ZZ_pE& x, const ZZ_p& a)
   { ::conv(x.rep, a); }

inline void ZZ_pEInfoT::conv(ZZ_pE& x, const ZZ& a)
   { ::conv(x.rep, a); }

inline ZZ_pE ZZ_pEInfoT::to_ZZ_pE(const ZZ_pX& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE ZZ_pEInfoT::to_ZZ_pE(long a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE ZZ_pEInfoT::to_ZZ_pE(const ZZ_p& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE ZZ_pEInfoT::to_ZZ_pE(const ZZ& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }


// ****** norm and trace

inline void ZZ_pEInfoT::trace(ZZ_p& x, const ZZ_pE& a)
   { TraceMod(x, a.rep, p); }
inline ZZ_p ZZ_pEInfoT::trace(const ZZ_pE& a)
   { return TraceMod(a.rep, p); }

inline void ZZ_pEInfoT::norm(ZZ_p& x, const ZZ_pE& a)
   { NormMod(x, a.rep, p); }
inline ZZ_p ZZ_pEInfoT::norm(const ZZ_pE& a)
   { return NormMod(a.rep, p); }


// ****** random numbers

inline void ZZ_pEInfoT::random(ZZ_pE& x)
// x = random element in ZZ_pE

   { ::random(x.rep, deg(p)); }

inline ZZ_pE ZZ_pEInfoT::random_ZZ_pE()
   { ZZ_pE x; random(x); NTL_OPT_RETURN(ZZ_pE, x); }


#endif
