/* -*- mode: c; style: linux -*- */

/* linbox/src/integer.inl
 * Copyright (C) 2001 B. David Saunders,
 *               2002 Jean-Guillaume Dumas,
 *               2002 Bradford Hovinen
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "linbox/integer.h"

namespace LinBox 
{

inline Integer::~Integer()
{
	mpz_clear((mpz_ptr)&gmp_rep);
}

inline Integer::Integer(const Integer &n) {
	mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
}

inline Integer& Integer::logcpy(const Integer &n)
{
	if (this == &n) return *this;
	mpz_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
	return *this;
}

// same as logcopy
inline Integer& Integer::operator = (const Integer &n) { return logcpy(n) ; }

//-----------------------------Integer(int n)
inline Integer::Integer(int n = 0) { mpz_init_set_si((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(uint n)
inline Integer::Integer(unsigned int n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(long n)
inline Integer::Integer(long n) { mpz_init_set_si((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(unsigned long n)
inline Integer::Integer(unsigned long n) { mpz_init_set_ui((mpz_ptr)&gmp_rep, n) ; }

//-----------------------------Integer(double)
inline Integer::Integer(double d) { mpz_init_set_d((mpz_ptr)&gmp_rep, d) ; }

//-------------------------------------------------inline comparaison operators
inline int operator != (const Integer& a , const Integer& b)
{ return compare(a,b) != 0; }

inline int operator != (int l, const Integer& n)
{ return n.operator != (l); }

inline int operator != (long l, const Integer& n)
{ return n.operator != (l); }

inline int operator == (const Integer& a, const Integer& b) 
{  return compare(a,b) == 0; }

inline int operator == (int l, const Integer& n)
{ return (! (n.operator != (l))); }

inline int operator == (long l, const Integer& n)
{ return (! (n.operator != (l))); }

inline int operator == (const Integer& n, int l)
{ return (! (n.operator != (l))); }

inline int operator == (const Integer& n, long l)
{ return (! (n.operator != (l))); }

inline int operator < (const Integer& a , const Integer& b)
{ return compare(a,b) < 0; }

inline int operator < (const int l, const Integer& n)
{ return n > l; }

inline int operator < (const long l, const Integer& n)
{ return n > l; }

inline int operator >  (const Integer& a , const Integer& b)
{ return compare(a,b) > 0; }

inline int operator > (int l, const Integer& n)
{ return n < l; }

inline int operator > (long l, const Integer& n)
{ return n < l; }

inline int operator <= (const Integer& a, const Integer& b)
{ return compare(a,b) <= 0; }

inline int operator <= (const Integer& n, int l)
{  return (! (n > l) ); }

inline int operator <= (const Integer& n, long l)
{  return (! (n > l) ); }

inline int operator <= (int l, const Integer& n)
{  return (! (n < l) );}

inline int operator <= (long l, const Integer& n)
{  return (! (n < l) );}

inline int operator >= (const Integer& a, const Integer& b)
{ return compare(a,b) >= 0; }

inline int operator >= (int l, const Integer& n)
{  return (! (n > l) );}

inline int operator >= (long l, const Integer& n)
{  return (! (n > l) );}

inline int operator >= (const Integer& n, int l)
{  return (! (n < l) );}

inline int operator >= (const Integer& n, long l)
{  return (! (n < l) );}


//----------------------------------arithmetic inline operators
inline Integer Integer::operator - () const
{
// JGD 18.06.1999
	Integer Res ;
	mpz_neg((mpz_ptr)&Res.gmp_rep, (mpz_ptr)&gmp_rep ); 
	return Res ;
}

// -- operator +
inline Integer operator + (const int l, const Integer& n) { return n + (long)l; }
inline Integer operator + (const unsigned int l, const Integer& n) { return n + (unsigned long)l; }
inline Integer operator + (const long l, const Integer& n) { return n + l; }
inline Integer operator + (const unsigned long l, const Integer& n) { return n + l; }
inline Integer operator + (const Integer& n, const int l) { return n + (long)l; }
inline Integer operator + (const Integer& n, const unsigned int l) { return n + (unsigned long)l; }

inline Integer& operator += (Integer& n, const int l) { return n += (long)l; }
inline Integer& operator += (Integer& n, const unsigned int l) { return n += (unsigned long)l; }

// -- operator -
inline Integer operator - (const int l, const Integer& n) { return -(n - (long)l); }
inline Integer operator - (const unsigned int l, const Integer& n) { return -(n - (unsigned long)l); }
inline Integer operator - (const long l, const Integer& n) { return -(n - l); }
inline Integer operator - (const unsigned long l, const Integer& n) { return -(n - l); }
inline Integer operator - (const Integer& n, const int l) { return n - (long)l; }
inline Integer operator - (const Integer& n, const unsigned int l) { return n - (unsigned long)l; }

inline Integer& operator -= (Integer& n, const int l) { return n -= (long)l; }
inline Integer& operator -= (Integer& n, const unsigned int l) { return n -= (unsigned long)l; }

// -- operator *
inline Integer operator * (const int l, const Integer& n) { return n * (long)l; }
inline Integer operator * (const unsigned int l, const Integer& n) { return n * (unsigned long)l; }
inline Integer operator * (const long l, const Integer& n) { return n * l; }
inline Integer operator * (const unsigned long l, const Integer& n) { return n * l; }
inline Integer operator * (const Integer& n, const int l) { return n * (long)l; }
inline Integer operator * (const Integer& n, const unsigned int l) { return n * (unsigned long)l; }

inline Integer& operator *= (Integer& n, const int l) { return n *= (long)l; }
inline Integer& operator *= (Integer& n, const unsigned int l) { return n *= (unsigned long)l; }

// -- operator /
inline Integer operator / (const int l, const Integer& n) { return Integer(l)/n; }
inline Integer operator / (const long l, const Integer& n) { return Integer(l)/n; }
inline Integer operator / (const Integer& n, const int l) { return n / (long)l; }
inline Integer operator / (const Integer& n, const unsigned int l)
{ return n / (unsigned long)l; }

inline Integer& operator /= (Integer& n, const int l)
{ if (l>=0) return n /= (unsigned long)l; else return  n = -(n / (unsigned long)-l); }
inline Integer& operator /= (Integer& n, const long l) { return n /= (unsigned long)l; }
inline Integer& operator /= (Integer& n, const unsigned int l) { return n /= (unsigned long)l; }

// -- operator %
inline Integer operator % (const int l, const Integer& n) { return Integer(l) % n; }
inline Integer operator % (const long l, const Integer& n) { return Integer(l) % n; }
inline Integer operator % (const Integer& n, const int l) { return n % (long)l; }
inline Integer operator % (const Integer& n, const unsigned int l) { return n % (unsigned long)l; }

inline Integer& operator %= (Integer& n, const int l) { return n %= (long)l; }
inline Integer& operator %= (Integer& n, const unsigned int l) { return n %= (unsigned long)l; }


//----------miscellaneous inline functions

inline int Integer::priv_sign() const { return mpz_sgn( (mpz_ptr)&gmp_rep ); }

inline int isone(const Integer& a) { return ! mpz_cmp_ui((mpz_ptr)&(a.gmp_rep), 1UL); }

inline int iszero(const Integer& a) { return ! mpz_cmp_ui((mpz_ptr)&(a.gmp_rep), 0UL); }

inline int iszero(const short int a) { return a ==0; }
inline int iszero(const int a) { return a ==0; }
inline int iszero(const long a) { return a ==0; }
inline int iszero(const unsigned short int a) { return a ==0; }
inline int iszero(const unsigned int a) { return a ==0; }
inline int iszero(const unsigned long a) { return a ==0; }

inline int sign(const Integer& a) { return a.priv_sign(); }

inline unsigned long length(const Integer& a) { return mpz_size( (mpz_ptr)&(a.gmp_rep) ) * sizeof(unsigned long); }

inline Integer abs(const Integer &n) { if (sign(n) >= 0) return n; return -n; }

inline size_t Integer::size() const { return  mpz_size( (mpz_ptr)&gmp_rep ) ; }

inline unsigned long Integer::operator[](size_t i) const
{
	if ( mpz_size( (mpz_ptr)&gmp_rep ) > i)
		return mpz_getlimbn( (mpz_ptr)&gmp_rep, i);
	else
		return 0;
}

inline Integer operator <<= (Integer& n, unsigned int l) {  return n = n << l; }
inline Integer operator <<= (Integer& n, unsigned long l) {  return n = n << l; } 
inline Integer operator >>= (Integer& n, unsigned int l) {  return n = n >> l; } 
inline Integer operator >>= (Integer& n, unsigned long l) {  return n = n >> l; }

//-------------------------------------------------inline >> & << operators
inline ostream& operator<< (ostream& o, const Integer& a) { return a.print(o); }
 
}
