/* -*- mode: c; style: linux -*- */

/* linbox/src/integer.h
 * Copyright (C) 2001 B. David Saunders,
 *               2002 Bradford Hovinen
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
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

#ifndef __INTEGER_H
#define __INTEGER_H

#include <vector>
#include <list>
#include <string>
#include <cstdlib>

extern "C" {
#    include "gmp.h"
}

#define SGN(n) (((n) == 0) ? 0 : (((n) < 0) ? -1 : 1))
#define ABS(n) (((n) < 0) ? -(n) : (n))

namespace LinBox
{
	class Integer {
	    public:
		typedef vector<mp_limb_t> vect_t;
		Integer (const vector<mp_limb_t> &vect_t);

		Integer (int n = 0);
		Integer (long n);
		Integer (unsigned int n);
		Integer (unsigned long n);
		Integer (double d);
		Integer (const char *s);
		Integer (const Integer &n);
		~Integer ();

		static const Integer zero;
		static const Integer one;

		// -- Assignment and copy operators
		Integer &operator = (const Integer &n);
		Integer &logcpy (const Integer &n);
		Integer &copy (const Integer &n);

		//------------------Equalities and inequalities between integers and longs
		int operator != (const int l) const;
		int operator != (const long l) const;
  
		friend int compare (const Integer &a, const Integer &b);
		friend int absCompare (const Integer &a, const Integer &b);

		int operator > (const int l) const;
		int operator > (const long l) const;
		int operator < (const int l) const;
		int operator < (const long l) const;

		//----------------Elementary arithmetic between Integers  &longs
		Integer &operator += (const Integer &n);  
		Integer &operator += (const unsigned long l);  
		Integer &operator += (const long l);  
		Integer  operator + (const Integer &n) const;  
		Integer  operator + (const unsigned long l) const;
		Integer  operator + (const long l) const;

		Integer &operator -= (const Integer &n);  
		Integer &operator -= (const unsigned long l);  
		Integer &operator -= (const long l);  
		Integer  operator - (const Integer &n) const;
		Integer  operator - (const unsigned long l) const;
		Integer  operator - (const long l) const;
		Integer  operator -() const;

		Integer &operator *= (const Integer &n);  
		Integer &operator *= (const unsigned long l);  
		Integer &operator *= (const long l);  
		Integer  operator * (const Integer &n) const;
		Integer  operator * (const unsigned long l) const;
		Integer  operator * (const long l) const;

		// -- Euclidian division of a/b: returns q or r such that
		// - a=b*q + r, with |r| < |b|, a*r >=0
		Integer &operator /= (const Integer &n);  
		Integer &operator /= (const unsigned long l);
		Integer &operator /= (const long l);
		Integer  operator /  (const Integer &n) const;
		Integer  operator /  (const unsigned long l) const;
		Integer  operator /  (const long l) const;

		Integer &operator %= (const Integer &n);  
		Integer &operator %= (const unsigned long l);
		Integer &operator %= (const long l);
		Integer  operator % (const Integer &n) const;
		long  operator % (const unsigned long l) const;
		long  operator % (const long l) const;

		// - Methods
		static Integer &addin (Integer &res, const Integer &n);  
		static Integer &addin (Integer &res, const long n);  
		static Integer &addin (Integer &res, const unsigned long n);  
		static Integer &add   (Integer &res, const Integer &n1, const Integer &n2);  
		static Integer &add   (Integer &res, const Integer &n1, const long n2);  
		static Integer &add   (Integer &res, const Integer &n1, const unsigned long n2);  

		static Integer &subin (Integer &res, const Integer &n);  
		static Integer &subin (Integer &res, const long n);  
		static Integer &subin (Integer &res, const unsigned long n);  
		static Integer &sub   (Integer &res, const Integer &n1, const Integer &n2);  
		static Integer &sub   (Integer &res, const Integer &n1, const long n2);  
		static Integer &sub   (Integer &res, const Integer &n1, const unsigned long n2);  
		static Integer &negin (Integer &res);  
		static Integer &neg   (Integer &res, const Integer &n);  

		static Integer &mulin (Integer &res, const Integer &n);  
		static Integer &mulin (Integer &res, const long n);  
		static Integer &mulin (Integer &res, const unsigned long n);  
		static Integer &mul   (Integer &res, const Integer &n1, const Integer &n2);  
		static Integer &mul   (Integer &res, const Integer &n1, const long n2);  
		static Integer &mul   (Integer &res, const Integer &n1, const unsigned long n2);  
		static Integer &axpy   (Integer &res, const Integer &a, const Integer &x, const Integer &y);  
		static Integer &axpyin   (Integer &res, const Integer &a, const Integer &x);  
		static Integer &axmy   (Integer &res, const Integer &a, const Integer &x, const Integer &y);  
		static Integer &axmyin   (Integer &res, const Integer &a, const Integer &x);  

		static Integer &divin (Integer &q, const Integer &n);  
		static Integer &divin (Integer &q, const long n);  
		static Integer &divin (Integer &q, const unsigned long n);  
		static Integer &div   (Integer &q, const Integer &n1, const Integer &n2);  
		static Integer &div   (Integer &q, const Integer &n1, const long n2);  
		static Integer &div   (Integer &q, const Integer &n1, const unsigned long n2);  
		static Integer &divexact  (Integer &q, const Integer &n1, const Integer &n2);  
		static Integer  divexact  (const Integer &n1, const Integer &n2);  

		static Integer &modin (Integer &r, const Integer &n);  
		static Integer &modin (Integer &r, const long n);  
		static Integer &modin (Integer &r, const unsigned long n);  
		static Integer &mod   (Integer &r, const Integer &n1, const Integer &n2);  
		static Integer &mod   (Integer &r, const Integer &n1, const long n2);  
		static Integer &mod   (Integer &r, const Integer &n1, const unsigned long n2);  

		// -- return q, the quotient
		static Integer &divmod   (Integer &q, Integer &r, const Integer &n1, const Integer &n2);  
		static Integer &divmod   (Integer &q, Integer &r, const Integer &n1, const long n2);  
		static Integer &divmod   (Integer &q, Integer &r, const Integer &n1, const unsigned long n2);  

  
		//------------------------------------- Arithmetic functions
		friend Integer gcd (const Integer &a, const Integer &b);
		friend Integer gcd (const Integer &a, const Integer &b, 
				    Integer &u, Integer &v);
		friend Integer &gcd (Integer &g, const Integer &a, const Integer &b);
		friend Integer &gcd (Integer &g, const Integer &a, const Integer &b, 
				     Integer &u, Integer &v);

		friend Integer pp (const Integer &P, const Integer &Q);


		// - return n^l 
		friend Integer pow (const Integer &n, const long l);
		friend Integer pow (const Integer &n, const unsigned long l);
		friend Integer pow (const Integer &n, const int l)
			{ return pow (n, (long) l); }
		friend Integer pow (const Integer &n, const unsigned int l)
			{ return pow (n, (unsigned long) l); }

		// - return n^e % m
		friend Integer powmod (const Integer &n, const unsigned long e, const Integer &m);
		friend Integer powmod (const Integer &n, const long e, const Integer &m);
		friend Integer powmod (const Integer &n, const unsigned int e, const Integer &m)
			{ return powmod (n, (unsigned long) e, m); }
		friend Integer powmod (const Integer &n, const int e, const Integer &m)
			{ return powmod (n, (long) e, m); }
		friend Integer powmod (const Integer &n, const Integer &e, const Integer &m);

		friend Integer fact ( unsigned long l);
  
		friend Integer sqrt (const Integer &p);
		friend Integer sqrt (const Integer &p, Integer &r);
		friend long logp (const Integer &a, const Integer &p) ;

		//-----------------------------------------Miscellaneous
		friend inline int sign   (const Integer &a);
		friend inline int iszero (const Integer &a);
		friend inline int isone  (const Integer &a);

		friend Integer abs (const Integer &n);

		friend int probab_prime (const Integer &p);
		friend int probab_prime (const Integer &p, int r);
		friend int jacobi (const Integer &u, const Integer &v) ;
		friend int legendre (const Integer &u, const Integer &v) ;


		Integer operator << (unsigned int l) const; // lshift
		Integer operator >> (unsigned int l) const; // rshift
		Integer operator << (unsigned long l) const; // lshift
		Integer operator >> (unsigned long l) const; // rshift

		// - return the size in byte
		friend inline unsigned long length (const Integer &a); 
		// - return the size in word.
		size_t size () const;
		// - return the i-th word of the integer. Word 0 is lowest word.
		unsigned long operator[](size_t i) const; 

		// -- Convert an Integer to a basic C++ type
		// -- Cast operators
		friend long Integer2long (const Integer &n);
		friend vect_t &Integer2vector  (vect_t &v, const Integer &n);
		friend double Integer2double (const Integer &n);
		friend string &Integer2string (string&, const Integer&, int base = 10);
		operator unsigned int () const ;
		operator int () const ;
		operator unsigned long () const ;
		operator long () const ;
		operator string () const ;
		operator float () const ;
		operator double () const ;
		operator vect_t () const ;


		//----------------------------------------------I/O

		friend istream &operator >> (istream &i, Integer &n);
		friend ostream &operator << (ostream &o, const Integer &n);
		friend ostream &absOutput (ostream &o, const Integer &n);

		ostream &print (ostream &o ) const;
  
	    protected:

		typedef MP_INT Rep;

		Rep gmp_rep;

		int priv_sign () const;

		// -- Creates a new Integer from a size sz and a array of unsigned long d 
		Integer (unsigned long *d, long size);

	}; //----------------------------------------------- End of Class Integer

	typedef int int_32;       // should guarantee 32 bits though
	typedef long long int_64; // should guarantee 64 bits though
}


#include "linbox/integer.inl"
#include "linbox/util/integer/integer-add.C"
#include "linbox/util/integer/integer-sub.C"
#include "linbox/util/integer/integer-mul.C"
#include "linbox/util/integer/integer-div.C"
#include "linbox/util/integer/integer-mod.C"
#include "linbox/util/integer/integer-pow.C"
#include "linbox/util/integer/integer-compare.C"
#include "linbox/util/integer/integer-io.C"
#include "linbox/util/integer/integer-misc.C"
#include "linbox/util/integer/integer-cstor.C"

#endif
