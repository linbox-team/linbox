
#ifndef NTL_zz_pE__H
#define NTL_zz_pE__H

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/vec_long.h>

#include <NTL/lzz_pX.h>


class zz_pE;

class zz_pEInfoT {
private:
   zz_pEInfoT();                       // disabled
   zz_pEInfoT(const zz_pEInfoT&);   // disabled
   void operator=(const zz_pEInfoT&);  // disabled
public:
   long ref_count;

   zz_pEInfoT(const zz_pX&);
   ~zz_pEInfoT() { }

   zz_pXModulus p;

   ZZ cardinality;

   void add(zz_pE& x, const zz_pE& a, const zz_pE& b);
   void sub(zz_pE& x, const zz_pE& a, const zz_pE& b);
   void negate(zz_pE& x, const zz_pE& a);

   void add(zz_pE& x, const zz_pE& a, long b);
   void add(zz_pE& x, const zz_pE& a, const zz_p& b);
   void add(zz_pE& x, long a, const zz_pE& b);
   void add(zz_pE& x, const zz_p& a, const zz_pE& b);

   void sub(zz_pE& x, const zz_pE& a, long b);
   void sub(zz_pE& x, const zz_pE& a, const zz_p& b);
   void sub(zz_pE& x, long a, const zz_pE& b);
   void sub(zz_pE& x, const zz_p& a, const zz_pE& b);

   void mul(zz_pE& x, const zz_pE& a, const zz_pE& b);
   void sqr(zz_pE& x, const zz_pE& a);
   zz_pE sqr(const zz_pE& a);

   void mul(zz_pE& x, const zz_pE& a, long b);
   void mul(zz_pE& x, const zz_pE& a, const zz_p& b);
   void mul(zz_pE& x, long a, const zz_pE& b);
   void mul(zz_pE& x, const zz_p& a, const zz_pE& b);

   void div(zz_pE& x, const zz_pE& a, const zz_pE& b);
   void div(zz_pE& x, const zz_pE& a, long b);
   void div(zz_pE& x, const zz_pE& a, const zz_p& b);
   void div(zz_pE& x, long a, const zz_pE& b);
   void div(zz_pE& x, const zz_p& a, const zz_pE& b);

   void inv(zz_pE& x, const zz_pE& a);
   zz_pE inv(const zz_pE& a);

   void power(zz_pE& x, const zz_pE& a, const ZZ& e);
   zz_pE power(const zz_pE& a, const ZZ& e);
   void power(zz_pE& x, const zz_pE& a, long e);
   zz_pE power(const zz_pE& a, long e);

   void conv(zz_pE& x, const zz_pX& a);
   void conv(zz_pE& x, long a);
   void conv(zz_pE& x, const zz_p& a);
   void conv(zz_pE& x, const ZZ& a);
   zz_pE to_zz_pE(const zz_pX& a);
   zz_pE to_zz_pE(long a);
   zz_pE to_zz_pE(const zz_p& a); 
   zz_pE to_zz_pE(const ZZ& a);

   void trace(zz_p& x, const zz_pE& a);
   zz_p trace(const zz_pE& a);

   void norm(zz_p& x, const zz_pE& a);
   zz_p norm(const zz_pE& a);

   void random(zz_pE& x);
   zz_pE random_zz_pE();

};

extern zz_pEInfoT *zz_pEInfo; // info for current modulus, initially null

#if (defined (_THREAD_SAFE)) || (defined (_REENTRANT))
#  if !defined (COARSE_LOCKS)

extern pthread_rwlock_t zz_pE_lock;

#  else
#    define zz_pE_lock field_lock;
#  endif
#  define NTL_LZZ_PE_THREADS_ENTER pthread_rwlock_rdlock (&zz_pE_lock);
#  define NTL_LZZ_PE_THREADS_LEAVE pthread_rwlock_unlock (&zz_pE_lock);
#else
#  define NTL_LZZ_PE_THREADS_ENTER
#  define NTL_LZZ_PE_THREADS_LEAVE
#endif




class zz_pEContext {
private:
zz_pEInfoT *ptr;

public:
void save();
void restore() const;

zz_pEContext() { ptr = 0; }
zz_pEContext(const zz_pX& p);

zz_pEContext(const zz_pEContext&); 


zz_pEContext& operator=(const zz_pEContext&); 


~zz_pEContext();


};


class zz_pEBak {
private:
long MustRestore;
zz_pEInfoT *ptr;

zz_pEBak(const zz_pEBak&); // disabled
void operator=(const zz_pEBak&); // disabled

public:
void save();
void restore();

zz_pEBak() { MustRestore = 0; ptr = 0; }

~zz_pEBak();


};



struct zz_pE_NoAlloc_type { zz_pE_NoAlloc_type() { } };
const zz_pE_NoAlloc_type zz_pE_NoAlloc = zz_pE_NoAlloc_type();



class zz_pE {


private:

friend class zz_pEInfoT;

zz_pX rep;

// static data


public:

static long DivCross() { return 16; }
static long ModCross() { return 8; }


// ****** constructors and assignment

zz_pE();

zz_pE(const zz_pE& a)  { rep.rep.SetMaxLength(zz_pE::degree()); rep = a.rep; }

zz_pE(zz_pE_NoAlloc_type) { }  // allocates no space

~zz_pE() { } 

zz_pE& operator=(const zz_pE& a) { rep = a.rep; return *this; }

zz_pE(zz_pE& x, INIT_TRANS_TYPE) : rep(x.rep, INIT_TRANS) { }


// read-only access to representation
friend const zz_pX& rep(const zz_pE& a) { return a.rep; }

// You can always access the representation directly...if you dare.
zz_pX& LoopHole() { return rep; }



static const zz_pXModulus& modulus() { return zz_pEInfo->p; }

static long degree() { return deg(zz_pEInfo->p); }

static const ZZ& cardinality() { return zz_pEInfo->cardinality; }

static const zz_pE& zero();

static long initialized() { return (zz_pEInfo != 0); }

static void init(const zz_pX&);

friend void clear(zz_pE& x)
// x = 0
   { clear(x.rep); }

friend void set(zz_pE& x)
// x = 1
   { set(x.rep); }

friend void swap(zz_pE& x, zz_pE& y)
// swap x and y

   { swap(x.rep, y.rep); }

// ****** addition

friend void add(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a + b

   { add(x.rep, a.rep, b.rep); }

friend void sub(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a - b

   { sub(x.rep, a.rep, b.rep); }


friend void negate(zz_pE& x, const zz_pE& a) 

   { negate(x.rep, a.rep); }


friend void add(zz_pE& x, const zz_pE& a, long b)
   { add(x.rep, a.rep, b); }

friend void add(zz_pE& x, const zz_pE& a, const zz_p& b)
   { add(x.rep, a.rep, b); }

friend void add(zz_pE& x, long a, const zz_pE& b)
   { add(x.rep, a, b.rep); }

friend void add(zz_pE& x, const zz_p& a, const zz_pE& b)
   { add(x.rep, a, b.rep); }





friend void sub(zz_pE& x, const zz_pE& a, long b)
   { sub(x.rep, a.rep, b); }

friend void sub(zz_pE& x, const zz_pE& a, const zz_p& b)
   { sub(x.rep, a.rep, b); }

friend void sub(zz_pE& x, long a, const zz_pE& b)
   { sub(x.rep, a, b.rep); }

friend void sub(zz_pE& x, const zz_p& a, const zz_pE& b)
   { sub(x.rep, a, b.rep); }





// ****** multiplication

friend void mul(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a*b

   { MulMod(x.rep, a.rep, b.rep, zz_pE::modulus()); }


friend void sqr(zz_pE& x, const zz_pE& a)
// x = a^2

   { SqrMod(x.rep, a.rep, zz_pE::modulus()); }

friend zz_pE sqr(const zz_pE& a)
   { zz_pE x; sqr(x, a); NTL_OPT_RETURN(zz_pE, x); }


friend void mul(zz_pE& x, const zz_pE& a, long b)
   { mul(x.rep, a.rep, b); }

friend void mul(zz_pE& x, const zz_pE& a, const zz_p& b)
   { mul(x.rep, a.rep, b); }

friend void mul(zz_pE& x, long a, const zz_pE& b)
   { mul(x.rep, a, b.rep); }

friend void mul(zz_pE& x, const zz_p& a, const zz_pE& b)
   { mul(x.rep, a, b.rep); }


// ****** division



friend void div(zz_pE& x, const zz_pE& a, const zz_pE& b) { zz_pEInfo->div (x, a, b); }
friend void div(zz_pE& x, const zz_pE& a, long b) { zz_pEInfo->div (x, a, b); }
friend void div(zz_pE& x, const zz_pE& a, const zz_p& b) { zz_pEInfo->div (x, a, b); }
friend void div(zz_pE& x, long a, const zz_pE& b) { zz_pEInfo->div (x, a, b); }
friend void div(zz_pE& x, const zz_p& a, const zz_pE& b) { zz_pEInfo->div (x, a, b); }

friend void inv(zz_pE& x, const zz_pE& a) { zz_pEInfo->inv (x, a); }
friend zz_pE inv(const zz_pE& a)
   { zz_pE x; inv(x, a); NTL_OPT_RETURN(zz_pE, x); }



// ****** exponentiation

friend void power(zz_pE& x, const zz_pE& a, const ZZ& e)
// x = a^e

   { PowerMod(x.rep, a.rep, e, zz_pE::modulus()); }

friend zz_pE power(const zz_pE& a, const ZZ& e)
   { zz_pE x; power(x, a, e); NTL_OPT_RETURN(zz_pE, x); }

friend void power(zz_pE& x, const zz_pE& a, long e)
   { _BUFFER ZZ ee; conv (ee, e); power(x, a, ee); }

friend zz_pE power(const zz_pE& a, long e)
   { zz_pE x; power(x, a, e); NTL_OPT_RETURN(zz_pE, x); }




// ****** conversion

friend void conv(zz_pE& x, const zz_pX& a)
   { rem(x.rep, a, zz_pE::modulus()); }

friend void conv(zz_pE& x, long a)
   { conv(x.rep, a); }

friend void conv(zz_pE& x, const zz_p& a)
   { conv(x.rep, a); }

friend void conv(zz_pE& x, const ZZ& a)
   { conv(x.rep, a); }

friend zz_pE to_zz_pE(const zz_pX& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

friend zz_pE to_zz_pE(long a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

friend zz_pE to_zz_pE(const zz_p& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

friend zz_pE to_zz_pE(const ZZ& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }



// ****** comparison

friend long IsZero(const zz_pE& a)
   { return IsZero(a.rep); }

friend long IsOne(const zz_pE& a)
   { return IsOne(a.rep); }

friend long operator==(const zz_pE& a, const zz_pE& b)
   { return a.rep == b.rep; }
friend long operator==(const zz_pE& a, long b)
   { return a.rep == b; }
friend long operator==(const zz_pE& a, const zz_p& b)
   { return a.rep == b; }
friend long operator==(long a, const zz_pE& b)
   { return a == b.rep; }
friend long operator==(const zz_p& a, const zz_pE& b)
   { return a == b.rep; }

friend long operator!=(const zz_pE& a, const zz_pE& b)
   { return !(a == b); }
friend long operator!=(const zz_pE& a, long b)
   { return !(a == b); }
friend long operator!=(const zz_pE& a, const zz_p& b)
   { return !(a == b); }
friend long operator!=(long a, const zz_pE& b)
   { return !(a == b); }
friend long operator!=(const zz_p& a, const zz_pE& b)
   { return !(a == b); }


// ****** norm and trace

friend void trace(zz_p& x, const zz_pE& a)
   { TraceMod(x, a.rep, zz_pE::modulus()); }
friend zz_p trace(const zz_pE& a)
   { return TraceMod(a.rep, zz_pE::modulus()); }

friend void norm(zz_p& x, const zz_pE& a)
   { NormMod(x, a.rep, zz_pE::modulus()); }
friend zz_p norm(const zz_pE& a)
   { return NormMod(a.rep, zz_pE::modulus()); }


// ****** random numbers

friend void random(zz_pE& x)
// x = random element in zz_pE

   { random(x.rep, zz_pE::degree()); }

friend zz_pE random_zz_pE()
   { zz_pE x; random(x); NTL_OPT_RETURN(zz_pE, x); }


// ****** input/output

friend ostream& operator<<(ostream& s, const zz_pE& a)
   { return s << a.rep; }
   
friend istream& operator>>(istream& s, zz_pE& x);


zz_pE& operator=(long a) { conv(*this, a); return *this; }
zz_pE& operator=(const zz_p& a) { conv(*this, a); return *this; }


};



inline zz_pE operator+(const zz_pE& a, const zz_pE& b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator+(const zz_pE& a, const zz_p& b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator+(const zz_pE& a, long b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator+(const zz_p& a, const zz_pE& b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator+(long a, const zz_pE& b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }


inline zz_pE operator-(const zz_pE& a, const zz_pE& b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(const zz_pE& a, const zz_p& b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(const zz_pE& a, long b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(const zz_p& a, const zz_pE& b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(long a, const zz_pE& b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(const zz_pE& a)
   { zz_pE x; negate(x, a); NTL_OPT_RETURN(zz_pE, x); } 


inline zz_pE& operator+=(zz_pE& x, const zz_pE& b)
   { add(x, x, b); return x; }

inline zz_pE& operator+=(zz_pE& x, const zz_p& b)
   { add(x, x, b); return x; }

inline zz_pE& operator+=(zz_pE& x, long b)
   { add(x, x, b); return x; }


inline zz_pE& operator-=(zz_pE& x, const zz_pE& b)
   { sub(x, x, b); return x; }

inline zz_pE& operator-=(zz_pE& x, const zz_p& b)
   { sub(x, x, b); return x; }

inline zz_pE& operator-=(zz_pE& x, long b)
   { sub(x, x, b); return x; }


inline zz_pE& operator++(zz_pE& x) { add(x, x, 1); return x; }

inline void operator++(zz_pE& x, int) { add(x, x, 1); }

inline zz_pE& operator--(zz_pE& x) { sub(x, x, 1); return x; }

inline void operator--(zz_pE& x, int) { sub(x, x, 1); }



inline zz_pE operator*(const zz_pE& a, const zz_pE& b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator*(const zz_pE& a, const zz_p& b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator*(const zz_pE& a, long b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator*(const zz_p& a, const zz_pE& b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator*(long a, const zz_pE& b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }


inline zz_pE& operator*=(zz_pE& x, const zz_pE& b)
   { mul(x, x, b); return x; }

inline zz_pE& operator*=(zz_pE& x, const zz_p& b)
   { mul(x, x, b); return x; }

inline zz_pE& operator*=(zz_pE& x, long b)
   { mul(x, x, b); return x; }




inline zz_pE operator/(const zz_pE& a, const zz_pE& b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator/(const zz_pE& a, const zz_p& b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator/(const zz_pE& a, long b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator/(const zz_p& a, const zz_pE& b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator/(long a, const zz_pE& b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }


inline zz_pE& operator/=(zz_pE& x, const zz_pE& b)
   { div(x, x, b); return x; }

inline zz_pE& operator/=(zz_pE& x, const zz_p& b)
   { div(x, x, b); return x; }

inline zz_pE& operator/=(zz_pE& x, long b)
   { div(x, x, b); return x; }





// Versions of the above operations that are methods of the class zz_pEInfoT




// ****** addition

inline void zz_pEInfoT::add(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a + b

   { ::add(x.rep, a.rep, b.rep); }

inline void zz_pEInfoT::sub(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a - b

   { ::sub(x.rep, a.rep, b.rep); }


inline void zz_pEInfoT::negate(zz_pE& x, const zz_pE& a) 

   { ::negate(x.rep, a.rep); }


inline void zz_pEInfoT::add(zz_pE& x, const zz_pE& a, long b)
   { ::add(x.rep, a.rep, b); }

inline void zz_pEInfoT::add(zz_pE& x, const zz_pE& a, const zz_p& b)
   { ::add(x.rep, a.rep, b); }

inline void zz_pEInfoT::add(zz_pE& x, long a, const zz_pE& b)
   { ::add(x.rep, a, b.rep); }

inline void zz_pEInfoT::add(zz_pE& x, const zz_p& a, const zz_pE& b)
   { ::add(x.rep, a, b.rep); }





inline void zz_pEInfoT::sub(zz_pE& x, const zz_pE& a, long b)
   { ::sub(x.rep, a.rep, b); }

inline void zz_pEInfoT::sub(zz_pE& x, const zz_pE& a, const zz_p& b)
   { ::sub(x.rep, a.rep, b); }

inline void zz_pEInfoT::sub(zz_pE& x, long a, const zz_pE& b)
   { ::sub(x.rep, a, b.rep); }

inline void zz_pEInfoT::sub(zz_pE& x, const zz_p& a, const zz_pE& b)
   { ::sub(x.rep, a, b.rep); }





// ****** multiplication

inline void zz_pEInfoT::mul(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a*b

   { MulMod(x.rep, a.rep, b.rep, p); }


inline void zz_pEInfoT::sqr(zz_pE& x, const zz_pE& a)
// x = a^2

   { SqrMod(x.rep, a.rep, p); }

inline zz_pE zz_pEInfoT::sqr(const zz_pE& a)
   { zz_pE x; sqr(x, a); NTL_OPT_RETURN(zz_pE, x); }


inline void zz_pEInfoT::mul(zz_pE& x, const zz_pE& a, long b)
   { ::mul(x.rep, a.rep, b); }

inline void zz_pEInfoT::mul(zz_pE& x, const zz_pE& a, const zz_p& b)
   { ::mul(x.rep, a.rep, b); }

inline void zz_pEInfoT::mul(zz_pE& x, long a, const zz_pE& b)
   { ::mul(x.rep, a, b.rep); }

inline void zz_pEInfoT::mul(zz_pE& x, const zz_p& a, const zz_pE& b)
   { ::mul(x.rep, a, b.rep); }


// ****** division

inline zz_pE zz_pEInfoT::inv(const zz_pE& a)
   { zz_pE x; inv(x, a); NTL_OPT_RETURN(zz_pE, x); }



// ****** exponentiation

inline void zz_pEInfoT::power(zz_pE& x, const zz_pE& a, const ZZ& e)
// x = a^e

   { PowerMod(x.rep, a.rep, e, p); }

inline zz_pE zz_pEInfoT::power(const zz_pE& a, const ZZ& e)
   { zz_pE x; power(x, a, e); NTL_OPT_RETURN(zz_pE, x); }

inline void zz_pEInfoT::power(zz_pE& x, const zz_pE& a, long e)
   { _BUFFER ZZ ee; ::conv (ee, e); power(x, a, ee); }

inline zz_pE zz_pEInfoT::power(const zz_pE& a, long e)
   { zz_pE x; power(x, a, e); NTL_OPT_RETURN(zz_pE, x); }




// ****** conversion

inline void zz_pEInfoT::conv(zz_pE& x, const zz_pX& a)
   { ::rem(x.rep, a, p); }

inline void zz_pEInfoT::conv(zz_pE& x, long a)
   { ::conv(x.rep, a); }

inline void zz_pEInfoT::conv(zz_pE& x, const zz_p& a)
   { ::conv(x.rep, a); }

inline void zz_pEInfoT::conv(zz_pE& x, const ZZ& a)
   { ::conv(x.rep, a); }

inline zz_pE zz_pEInfoT::to_zz_pE(const zz_pX& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE zz_pEInfoT::to_zz_pE(long a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE zz_pEInfoT::to_zz_pE(const zz_p& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE zz_pEInfoT::to_zz_pE(const ZZ& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }


// ****** norm and trace

inline void zz_pEInfoT::trace(zz_p& x, const zz_pE& a)
   { TraceMod(x, a.rep, p); }
inline zz_p zz_pEInfoT::trace(const zz_pE& a)
   { return TraceMod(a.rep, p); }

inline void zz_pEInfoT::norm(zz_p& x, const zz_pE& a)
   { NormMod(x, a.rep, p); }
inline zz_p zz_pEInfoT::norm(const zz_pE& a)
   { return NormMod(a.rep, p); }


// ****** random numbers

inline void zz_pEInfoT::random(zz_pE& x)
// x = random element in zz_pE

   { ::random(x.rep, deg (p)); }

inline zz_pE zz_pEInfoT::random_zz_pE()
   { zz_pE x; random(x); NTL_OPT_RETURN(zz_pE, x); }

#endif
