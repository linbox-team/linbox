/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/util/integer/integer-misc.C
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

#include <iostream>
#include <list>

#include "linbox/integer.h"

namespace LinBox 
{

//-------------------------------------------fact (unsigned long l)
integer fact ( unsigned long l) 
{
	integer Res ;
	mpz_fac_ui( (mpz_ptr)&(Res.gmp_rep), l ) ;
	return Res ;
}

//-------------------------------------------square root
integer sqrt(const integer &a)
{
	integer q;
	mpz_sqrt( (mpz_ptr)&(q.gmp_rep),
		  (mpz_ptr)&(a.gmp_rep)) ;
	return q;
}

integer sqrt(const integer &a, integer& r)
{
	integer q;
	mpz_sqrtrem( (mpz_ptr)&(q.gmp_rep),
		     (mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(a.gmp_rep)) ;
	return q;
}


// base p logarithm of a
long logp(const integer& a, const integer& p) {
	list< integer > pows;
	integer puiss = p, sq;
	do {
		pows.push_back( puiss );
	} while ( (puiss *= puiss) <= a );
	puiss = pows.back(); pows.pop_back();
	long res = (1 << pows.size());
	while (! pows.empty() ) {
		if ((sq = puiss * pows.back()) <= a) {
			puiss = sq;
			pows.pop_back();
			res += (1 << pows.size());
		} else
			pows.pop_back();
	}
	return res;
}
    

//------------------------------------------GMP isprime
//     If this function returns 0, OP is definitely not prime.  If it
//     returns 1, then OP is `probably' prime.  The probability of a
//     false positive is (1/4)^r.  A reasonable value of r is 25.

int probab_prime(const integer &p)
{
	return mpz_probab_prime_p ((mpz_ptr)&(p.gmp_rep),1) ;
}

int probab_prime(const integer &p, int r)
{
	return mpz_probab_prime_p ((mpz_ptr)&(p.gmp_rep),r) ;
}

// ==========================================================================
// Computes and returns the Jacobi and Legendre symbols (u/v) of the integers u and v.  
// The algorithm used is Gmp's.
int jacobi(const integer& u, const integer& v)
{
	return mpz_jacobi ((mpz_ptr)&(u.gmp_rep),(mpz_ptr)&(v.gmp_rep)) ;
}

int legendre(const integer& u, const integer& v)
{
	return mpz_legendre ((mpz_ptr)&(u.gmp_rep),(mpz_ptr)&(v.gmp_rep)) ;
}



//--------------------------------------------integer::operator <<   // N O T   I M P L E M E N T E D 
integer integer::operator << (unsigned int l) const 
{ return *this; }

integer integer::operator << (unsigned long l) const 
{ return *this; }


//--------------------------------------------integer::operator >>   // N O T   I M P L E M E N T E D 
integer integer::operator >> (unsigned int l) const
{ return *this; }

integer integer::operator >> (unsigned long l) const
{ return *this; }

//------------------------------------------- convert method
long integer2long  ( const integer& n)
{
	return mpz_get_si ( (mpz_srcptr)&n.gmp_rep);
}
double integer2double( const integer& n)
{
	return mpz_get_d( (mpz_srcptr)&n.gmp_rep);
}

integer::operator long() const {
	return mpz_get_si ( (mpz_srcptr)&gmp_rep);
}
integer::operator double() const {
	return mpz_get_d ( (mpz_srcptr)&gmp_rep);
}

integer::operator int() const {
        return mpz_get_si ( (mpz_srcptr)&gmp_rep);
}

integer::operator unsigned int() const {
        return mpz_get_ui ( (mpz_srcptr)&gmp_rep);
}

integer::operator unsigned long() const {
        return mpz_get_ui ( (mpz_srcptr)&gmp_rep);
}

integer::operator float() const {
        return (float)mpz_get_d ( (mpz_srcptr)&gmp_rep);
}

 
}

