/* -*- mode: c; style: linux -*- */

/* linbox/integer.h
 *
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
 *                        JGD Random functions back.                          
 *                        (2002/02/12 16:05:24) 
 *
 */

#ifndef __INTEGER_H
#define __INTEGER_H

#include <vector>
#include <list>
#include <string>
#include <cstdlib>


extern "C" {
#include "gmp.h"
}

#define SGN(n) (((n) == 0) ? 0 : (((n) < 0) ? -1 : 1))
#define ABS(n) (((n) < 0) ? -(n) : (n))


namespace LinBox
{
        using namespace std;

	class integer {
	    public:
		typedef vector<mp_limb_t> vect_t;
		integer (const vector<mp_limb_t> &vect_t);

		integer (int n = 0);
		integer (long n);
		integer (unsigned int n);
		integer (unsigned long n);
		integer (double d);
		integer (const char *s);
		integer (const integer &n);
		~integer ();

		static const integer zero;
		static const integer one;

		// -- Assignment and copy operators
		integer &operator = (const integer &n);
		integer &logcpy (const integer &n);
		integer &copy (const integer &n);

		//------------------Equalities and inequalities between integers and longs
		int operator != (const int l) const;
		int operator != (const long l) const;
  
		friend int compare (const integer &a, const integer &b);
		friend int absCompare (const integer &a, const integer &b);

		int operator > (const int l) const;
		int operator > (const long l) const;
		int operator < (const int l) const;
		int operator < (const long l) const;

		//----------------elementary arithmetic between integers  &longs
		integer &operator += (const integer &n);  
		integer &operator += (const unsigned long l);  
		integer &operator += (const long l);  
		integer  operator + (const integer &n) const;  
		integer  operator + (const unsigned long l) const;
		integer  operator + (const long l) const;

		integer &operator -= (const integer &n);  
		integer &operator -= (const unsigned long l);  
		integer &operator -= (const long l);  
		integer  operator - (const integer &n) const;
		integer  operator - (const unsigned long l) const;
		integer  operator - (const long l) const;
		integer  operator -() const;

		integer &operator *= (const integer &n);  
		integer &operator *= (const unsigned long l);  
		integer &operator *= (const long l);  
		integer  operator * (const integer &n) const;
		integer  operator * (const unsigned long l) const;
		integer  operator * (const long l) const;

		// -- Euclidian division of a/b: returns q or r such that
		// - a=b*q + r, with |r| < |b|, a*r >=0
		integer &operator /= (const integer &n);  
		integer &operator /= (const unsigned long l);
		integer &operator /= (const long l);
		integer  operator /  (const integer &n) const;
		integer  operator /  (const unsigned long l) const;
		integer  operator /  (const long l) const;

		integer &operator %= (const integer &n);  
		integer &operator %= (const unsigned long l);
		integer &operator %= (const long l);
		integer  operator % (const integer &n) const;
		long  operator % (const unsigned long l) const;
		long  operator % (const long l) const;

		// - Methods
		static integer &addin (integer &res, const integer &n);  
		static integer &addin (integer &res, const long n);  
		static integer &addin (integer &res, const unsigned long n);  
		static integer &add   (integer &res, const integer &n1, const integer &n2);  
		static integer &add   (integer &res, const integer &n1, const long n2);  
		static integer &add   (integer &res, const integer &n1, const unsigned long n2);  

		static integer &subin (integer &res, const integer &n);  
		static integer &subin (integer &res, const long n);  
		static integer &subin (integer &res, const unsigned long n);  
		static integer &sub   (integer &res, const integer &n1, const integer &n2);  
		static integer &sub   (integer &res, const integer &n1, const long n2);  
		static integer &sub   (integer &res, const integer &n1, const unsigned long n2);  
		static integer &negin (integer &res);  
		static integer &neg   (integer &res, const integer &n);  

		static integer &mulin (integer &res, const integer &n);  
		static integer &mulin (integer &res, const long n);  
		static integer &mulin (integer &res, const unsigned long n);  
		static integer &mul   (integer &res, const integer &n1, const integer &n2);  
		static integer &mul   (integer &res, const integer &n1, const long n2);  
		static integer &mul   (integer &res, const integer &n1, const unsigned long n2);  
		static integer &axpy   (integer &res, const integer &a, const integer &x, const integer &y);  
		static integer &axpyin   (integer &res, const integer &a, const integer &x);  
		static integer &axmy   (integer &res, const integer &a, const integer &x, const integer &y);  
		static integer &axmyin   (integer &res, const integer &a, const integer &x);  

		static integer &divin (integer &q, const integer &n);  
		static integer &divin (integer &q, const long n);  
		static integer &divin (integer &q, const unsigned long n);  
		static integer &div   (integer &q, const integer &n1, const integer &n2);  
		static integer &div   (integer &q, const integer &n1, const long n2);  
		static integer &div   (integer &q, const integer &n1, const unsigned long n2);  
		static integer &divexact  (integer &q, const integer &n1, const integer &n2);  
		static integer  divexact  (const integer &n1, const integer &n2);  

		static integer &modin (integer &r, const integer &n);  
		static integer &modin (integer &r, const long n);  
		static integer &modin (integer &r, const unsigned long n);  
		static integer &mod   (integer &r, const integer &n1, const integer &n2);  
		static integer &mod   (integer &r, const integer &n1, const long n2);  
		static integer &mod   (integer &r, const integer &n1, const unsigned long n2);  

		// -- return q, the quotient
		static integer &divmod   (integer &q, integer &r, const integer &n1, const integer &n2);  
		static integer &divmod   (integer &q, integer &r, const integer &n1, const long n2);  
		static integer &divmod   (integer &q, integer &r, const integer &n1, const unsigned long n2);  

  
		//------------------------------------- Arithmetic functions
		friend integer gcd (const integer &a, const integer &b);
		friend integer gcd (const integer &a, const integer &b, 
				    integer &u, integer &v);
		friend integer &gcd (integer &g, const integer &a, const integer &b);
		friend integer &gcd (integer &g, const integer &a, const integer &b, 
				     integer &u, integer &v);

		friend integer pp (const integer &P, const integer &Q);


		// - return n^l 
		friend integer pow (const integer &n, const long l);
		friend integer pow (const integer &n, const unsigned long l);
		friend integer pow (const integer &n, const int l)
			{ return pow (n, (long) l); }
		friend integer pow (const integer &n, const unsigned int l)
			{ return pow (n, (unsigned long) l); }

		// - return n^e % m
		friend integer powmod (const integer &n, const unsigned long e, const integer &m);
		friend integer powmod (const integer &n, const long e, const integer &m);
		friend integer powmod (const integer &n, const unsigned int e, const integer &m)
			{ return powmod (n, (unsigned long) e, m); }
		friend integer powmod (const integer &n, const int e, const integer &m)
			{ return powmod (n, (long) e, m); }
		friend integer powmod (const integer &n, const integer &e, const integer &m);

		friend integer fact ( unsigned long l);
  
		friend integer sqrt (const integer &p);
		friend integer sqrt (const integer &p, integer &r);
		friend long logp (const integer &a, const integer &p) ;

		//-----------------------------------------Miscellaneous
		friend inline int sign   (const integer &a);
		friend inline int iszero (const integer &a);
		friend inline int isone  (const integer &a);

		friend integer abs (const integer &n);

		friend int probab_prime (const integer &p);
		friend int probab_prime (const integer &p, int r);
		friend int jacobi (const integer &u, const integer &v) ;
		friend int legendre (const integer &u, const integer &v) ;


		integer operator << (unsigned int l) const; // lshift
		integer operator >> (unsigned int l) const; // rshift
		integer operator << (unsigned long l) const; // lshift
		integer operator >> (unsigned long l) const; // rshift

		// - return the size in byte
		friend inline unsigned long length (const integer &a); 
		// - return the size in word.
		size_t size () const;
		// - return the i-th word of the integer. Word 0 is lowest word.
		unsigned long operator[](size_t i) const; 

		// -- Convert an integer to a basic C++ type
		// -- Cast operators
		friend long integer2long (const integer &n);
		friend vect_t &integer2vector  (vect_t &v, const integer &n);
		friend double integer2double (const integer &n);
		friend string &integer2string (string&, const integer&, int base = 10);
		operator unsigned int () const ;
		operator int () const ;
		operator unsigned long () const ;
		operator long () const ;
		operator string () const ;
		operator float () const ;
		operator double () const ;
		operator vect_t () const ;

		//------------------Random Iterators
		// -- return a random number with sz machine word.
		// -- To be improved.
		static integer  random(int sz=1 );
		static integer  nonzerorandom(int sz=1 );
		static integer& random(integer& r, const integer& size );
		static integer& nonzerorandom(integer& r, 
                                             const integer& size );
		static integer& random(integer& r, long size =1 );
		static integer& nonzerorandom(integer& r, long size =1 );


		//----------------------------------------------I/O

		friend istream &operator >> (istream &i, integer &n);
		friend ostream &operator << (ostream &o, const integer &n);
		friend ostream &absOutput (ostream &o, const integer &n);

		ostream &print (ostream &o ) const;
  
	    protected:

		typedef MP_INT Rep;

		Rep gmp_rep;

		int priv_sign () const;

		// -- Creates a new integer from a size sz and a array of unsigned long d 
		integer (unsigned long *d, long size);

	}; //----------------------------------------------- End of Class integer

	typedef int int_32;       // should guarantee 32 bits though
	typedef long long int_64; // should guarantee 64 bits though
}


#include "linbox/integer.inl"

#endif
