/* -*- mode: c; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/integer.inl
 *
 * Copyright(c)'94-97 by Givaro Team
 * Copyright(c)'2000-2002 by LinBox Team
 * see the copyright file.
 * Created by M. Samama, T. Gautier
 *
 * Modified Jean-Guillaume.Dumas <Jean-Guillaume.Dumas@imag.fr>
 *          B. David Saunders <saunders@cis.udel.edu>,
 *          Bradford Hovinen <hovinen@cis.udel.edu>
 *          Gilles Villard <Gilles.Villard@ens-lyon.fr>
 *                        Fri Apr  5 16:55:46 EST 2002
 *                        JGD Random functions back.
 *                        (2002/02/12 16:05:24) 
 *
 */

#include "linbox/integer.h"

namespace LinBox 
{

inline integer::~integer()
{
	mpz_clear((mpz_ptr)&gmp_rep);
}

inline integer::integer(const integer &n) {
	mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
}

inline integer& integer::logcpy(const integer &n)
{
	if (this == &n) return *this;
	mpz_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
	return *this;
}

// same as logcopy
inline integer& integer::operator = (const integer &n) { return logcpy(n) ; }

//-----------------------------integer(int n)
inline integer::integer(int n = 0) { mpz_init_set_si((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------integer(uint n)
inline integer::integer(unsigned int n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------integer(long n)
inline integer::integer(long n) { mpz_init_set_si((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------integer(unsigned long n)
inline integer::integer(unsigned long n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------integer(double)
inline integer::integer(double d) { mpz_init_set_d((mpz_ptr)&gmp_rep, d) ; }

//-------------------------------------------------inline comparaison operators
inline int operator != (const integer& a , const integer& b)
{ return compare(a,b) != 0; }

inline int operator != (int l, const integer& n)
{ return n.operator != (l); }

inline int operator != (long l, const integer& n)
{ return n.operator != (l); }

inline int operator == (const integer& a, const integer& b) 
{  return compare(a,b) == 0; }

inline int operator == (int l, const integer& n)
{ return (! (n.operator != (l))); }

inline int operator == (long l, const integer& n)
{ return (! (n.operator != (l))); }

inline int operator == (const integer& n, int l)
{ return (! (n.operator != (l))); }

inline int operator == (const integer& n, long l)
{ return (! (n.operator != (l))); }

inline int operator < (const integer& a , const integer& b)
{ return compare(a,b) < 0; }

inline int operator < (const int l, const integer& n)
{ return n > l; }

inline int operator < (const long l, const integer& n)
{ return n > l; }

inline int operator >  (const integer& a , const integer& b)
{ return compare(a,b) > 0; }

inline int operator > (int l, const integer& n)
{ return n < l; }

inline int operator > (long l, const integer& n)
{ return n < l; }

inline int operator <= (const integer& a, const integer& b)
{ return compare(a,b) <= 0; }

inline int operator <= (const integer& n, int l)
{  return (! (n > l) ); }

inline int operator <= (const integer& n, long l)
{  return (! (n > l) ); }

inline int operator <= (int l, const integer& n)
{  return (! (n < l) );}

inline int operator <= (long l, const integer& n)
{  return (! (n < l) );}

inline int operator >= (const integer& a, const integer& b)
{ return compare(a,b) >= 0; }

inline int operator >= (int l, const integer& n)
{  return (! (n > l) );}

inline int operator >= (long l, const integer& n)
{  return (! (n > l) );}

inline int operator >= (const integer& n, int l)
{  return (! (n < l) );}

inline int operator >= (const integer& n, long l)
{  return (! (n < l) );}


//----------------------------------arithmetic inline operators
inline integer integer::operator - () const
{
// JGD 18.06.1999
	integer Res ;
	mpz_neg((mpz_ptr)&Res.gmp_rep, (mpz_ptr)&gmp_rep ); 
	return Res ;
}

// -- operator +
inline integer operator + (const int l, const integer& n) { return n + (long)l; }
inline integer operator + (const unsigned int l, const integer& n) { return n + (unsigned long)l; }
inline integer operator + (const long l, const integer& n) { return n + l; }
inline integer operator + (const unsigned long l, const integer& n) { return n + l; }
inline integer operator + (const integer& n, const int l) { return n + (long)l; }
inline integer operator + (const integer& n, const unsigned int l) { return n + (unsigned long)l; }

inline integer& operator += (integer& n, const int l) { return n += (long)l; }
inline integer& operator += (integer& n, const unsigned int l) { return n += (unsigned long)l; }

// -- operator -
inline integer operator - (const int l, const integer& n) { return -(n - (long)l); }
inline integer operator - (const unsigned int l, const integer& n) { return -(n - (unsigned long)l); }
inline integer operator - (const long l, const integer& n) { return -(n - l); }
inline integer operator - (const unsigned long l, const integer& n) { return -(n - l); }
inline integer operator - (const integer& n, const int l) { return n - (long)l; }
inline integer operator - (const integer& n, const unsigned int l) { return n - (unsigned long)l; }

inline integer& operator -= (integer& n, const int l) { return n -= (long)l; }
inline integer& operator -= (integer& n, const unsigned int l) { return n -= (unsigned long)l; }

// -- operator *
inline integer operator * (const int l, const integer& n) { return n * (long)l; }
inline integer operator * (const unsigned int l, const integer& n) { return n * (unsigned long)l; }
inline integer operator * (const long l, const integer& n) { return n * l; }
inline integer operator * (const unsigned long l, const integer& n) { return n * l; }
inline integer operator * (const integer& n, const int l) { return n * (long)l; }
inline integer operator * (const integer& n, const unsigned int l) { return n * (unsigned long)l; }

inline integer& operator *= (integer& n, const int l) { return n *= (long)l; }
inline integer& operator *= (integer& n, const unsigned int l) { return n *= (unsigned long)l; }

// -- operator /
inline integer operator / (const int l, const integer& n) { return integer(l)/n; }
inline integer operator / (const long l, const integer& n) { return integer(l)/n; }
inline integer operator / (const integer& n, const int l) { return n / (long)l; }
inline integer operator / (const integer& n, const unsigned int l)
{ return n / (unsigned long)l; }

inline integer& operator /= (integer& n, const int l)
{ if (l>=0) return n /= (unsigned long)l; else return  n = -(n / (unsigned long)-l); }
inline integer& operator /= (integer& n, const long l) { return n /= (unsigned long)l; }
inline integer& operator /= (integer& n, const unsigned int l) { return n /= (unsigned long)l; }

// -- operator %
inline integer operator % (const int l, const integer& n) { return integer(l) % n; }
inline integer operator % (const long l, const integer& n) { return integer(l) % n; }
inline integer operator % (const integer& n, const int l) { return n % (long)l; }
inline integer operator % (const integer& n, const unsigned int l) { return n % (unsigned long)l; }

inline integer& operator %= (integer& n, const int l) { return n %= (long)l; }
inline integer& operator %= (integer& n, const unsigned int l) { return n %= (unsigned long)l; }


//----------miscellaneous inline functions

inline int integer::priv_sign() const { return mpz_sgn( (mpz_ptr)&gmp_rep ); }

inline int isone(const integer& a) { return ! mpz_cmp_ui((mpz_ptr)&(a.gmp_rep), 1UL); }

inline int iszero(const integer& a) { return ! mpz_cmp_ui((mpz_ptr)&(a.gmp_rep), 0UL); }

inline int iszero(const short int a) { return a ==0; }
inline int iszero(const int a) { return a ==0; }
inline int iszero(const long a) { return a ==0; }
inline int iszero(const unsigned short int a) { return a ==0; }
inline int iszero(const unsigned int a) { return a ==0; }
inline int iszero(const unsigned long a) { return a ==0; }

inline int sign(const integer& a) { return a.priv_sign(); }

inline unsigned long length(const integer& a) { return mpz_size( (mpz_ptr)&(a.gmp_rep) ) * sizeof(unsigned long); }

inline integer abs(const integer &n) { if (sign(n) >= 0) return n; return -n; }

inline size_t integer::size() const { return  mpz_size( (mpz_ptr)&gmp_rep ) ; }

inline unsigned long integer::operator[](size_t i) const
{
	if ( mpz_size( (mpz_ptr)&gmp_rep ) > i)
		return mpz_getlimbn( (mpz_ptr)&gmp_rep, i);
	else
		return 0;
}

inline integer operator <<= (integer& n, unsigned int l) {  return n = n << l; }
inline integer operator <<= (integer& n, unsigned long l) {  return n = n << l; } 
inline integer operator >>= (integer& n, unsigned int l) {  return n = n >> l; } 
inline integer operator >>= (integer& n, unsigned long l) {  return n = n >> l; }

//-------------------------------------------------inline >> & << operators
inline ostream& operator<< (ostream& o, const integer& a) { return a.print(o); }
 

//----------------------- Random integers ----------

inline integer integer::random(int sz = 1)
{
  integer res;
  mpz_random((mpz_ptr) &(res.gmp_rep), sz);
  return res;
}

inline integer integer::nonzerorandom(int sz) {
    integer r;
    while(iszero(r  = random(sz) )) {};
    return r;
}

inline integer& integer::random (integer& r, const integer& similar)
{
     mpz_random((mpz_ptr) &(r.gmp_rep), mpz_size( (mpz_ptr)&(similar.gmp_rep) )
);
     mpz_tdiv_r( (mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(similar.gmp_rep) );
     return r;
};

inline integer& integer::nonzerorandom (integer& r, const integer& size) {
    while (iszero(r = random(r,size))) {};
    return r;
}


inline integer& integer::random (integer& r, long size = 1)
{
    mpz_random((mpz_ptr) &(r.gmp_rep), size);
    return r;
};


inline integer& integer::nonzerorandom (integer& r, long size = 1)
{    while (iszero(r = random(r,size))) {};
    return r;
}


}
