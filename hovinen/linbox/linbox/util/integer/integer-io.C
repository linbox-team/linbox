/* -*- mode: c; style: linux -*- */

/* linbox/linbox/util/integer/integer-io.C
 * Copyright (C) Givaro Team
 *
 * Written by M. Samama, T. Gautier
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

#include <iostream.h>
#include <stdlib.h>

#include "linbox/integer.h"

namespace LinBox
{

// Sortie nonsignee : 321321 meme si n = -321321, par exemple 
ostream& absOutput(ostream &o, const Integer&n)
{
	int base = 10;
  
	unsigned long strSize = mpz_sizeinbase((mpz_ptr)&(n.gmp_rep), base) + 2;
	char *str = ::new char[strSize];
	mpz_get_str(str, base, (mpz_ptr)&(n.gmp_rep));
	if (sign(n) < 0) {
		char *str1 = &(str[1]) ;
		o << str1;
	}
	else o << str;
	delete [] str ;
	return o;
} 

// Sortie signee : +321321 ou -321321, par exemple 
ostream& Integer::print(ostream &o) const
{
	// -- bug ?
	int base = 10;
	unsigned long strSize = mpz_sizeinbase((mpz_ptr)&(gmp_rep), base) + 2;
	char *str = new char[strSize];
	mpz_get_str(str, base, (mpz_ptr)&(gmp_rep));
// JGD 08.11.1999 : temporaire
//   if (sign(*this) > 0) o << '+' ;
	o << str;
	delete [] str ;
	return o;
} 

string& Integer2string(string& s, const Integer& n, int base = 10) {
	unsigned long strSize = mpz_sizeinbase((mpz_ptr)&(n.gmp_rep), base) + 2;
	char *str = new char[strSize + 2];
	mpz_get_str(str, base, (mpz_ptr)&(n.gmp_rep));
	s = string(str);
//    delete [] str ;
	return s;
}
Integer::operator string () const {
	string s;
	return Integer2string(s,*this);
}


Integer::Integer(const vector<mp_limb_t>& v) {
 	size_t s = v.size();
	if (s) {
	 	mpz_init_set_ui((mpz_ptr)&gmp_rep, v[0]);
		Integer base(256), prod;
		prod = base = pow(base, (unsigned long)sizeof(mp_limb_t) );

		vector<mp_limb_t>::const_iterator vi = v.begin();
		for(++vi;vi != v.end();++vi) { 
	 		*this += ( prod * (*vi) );
			prod *= base;
		}
	} else
	 	mpz_init( (mpz_ptr)&gmp_rep );

}

vector<mp_limb_t>& Integer2vector(vector<mp_limb_t>& v, const Integer& n) {
	size_t s = mpz_size( (mpz_ptr)&n.gmp_rep );
	v.resize(s); 
	vector<mp_limb_t>::iterator vi = v.begin();
	for(mp_size_t i = 0;vi != v.end();++vi, ++i) *vi = mpz_getlimbn( (mpz_ptr)&n.gmp_rep ,i);
}

Integer::operator vector<mp_limb_t> () const {
	vector<mp_limb_t> v;
	return Integer2vector(v,*this);
}

// Entree au format de la sortie
istream& operator>> (istream& in, Integer& a)
{
	static long base[10] = {
		10,
		100,
		1000,
		10000,
		100000,
		1000000,
		10000000,
		100000000,
		1000000000
	} ;
	if (!in.good())
	{
//		throw GivError("*** error: corruped stream, in operator>> for Integer") ;
	}
	if (!in) return in ;
	// eat white
	in >> ws  ;

	// Base : 10^9, we read by packet of length 9
	// the char.
	char Tmp[10] ;
	int counter = 0 ;

	// Set the returned integer
	a = 0L ;
	char ch ;
	int sign = 1 ;

	// find a sign:
	in.get(ch) ;
	if ((ch != '+') && (ch != '-') && !((ch >= '0') && (ch <= '9')))
	{
		cerr << "Bad integer format: found: "<< ch ;
		cerr << ", in place of '+' '-' or a digit"<< endl ;
		return in ;
	}
	switch (ch) {
	case '+' : break ;
	case '-' : sign = -1 ; break ;
	default  : in.putback(ch) ; break ;
	}
	// eat white
	in >> ws  ;

	int noend = 1 ;
	while (noend)
	{
		counter = 0 ;

		// Read 9 digits or less
		while ((noend) && (counter < 9)) {
			in.get(ch) ;
			if (in.eof()) { noend = 0 ; }
			else if ((ch >= '0') && (ch <= '9')) Tmp[counter++] = ch ;
			else { noend = 0 ;  in.putback(ch) ; }
		}
		if (counter >0) {
			long l ;
			Tmp[counter] = '\0' ; // terminate the string
			l = atol(Tmp) ;
			a = a * base[counter-1] + l ;
		}
	}
	if (sign == -1) a = -a ;  
	return in ;
}
 
}
