/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/util/integer/integer-mod.C
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

#include "linbox/integer.h"
#include "linbox/util/error.h"

namespace LinBox 
{

//-------------------------------------------------- operator /
integer& integer::modin(integer& res, const integer& n) 
{
	if (iszero(n)) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(res)) return res;
	mpz_tdiv_r( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
	return res;
}
integer& integer::modin(integer& res, const long n) 
{
	if (n ==0) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(res)) return res;
	int sgn = SGN(n); 
	mpz_tdiv_r_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, ABS(n));
	if (sgn <0) return res = -res;
	return res;
}
integer& integer::modin(integer& res, const unsigned long n) 
{
	if (n ==0) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(res)) return res;
	mpz_tdiv_r_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
	return res;
}

integer& integer::mod(integer& res, const integer& n1, const integer& n2)
{
	if (iszero(n1)) return res = integer::zero;
	if (iszero(n2)) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	mpz_tdiv_r( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
	return res;
}
integer& integer::mod(integer& res, const integer& n1, const long n2)
{
	if (iszero(n1)) return res = integer::zero;
	if (iszero(n2)) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	int sgn = SGN(n2); 
	mpz_tdiv_r_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, ABS(n2));
	if (sgn <0) return res = - res;
	return res;
}
integer& integer::mod(integer& res, const integer& n1, const unsigned long n2)
{
	if (iszero(n1)) return res = integer::zero;
	if (iszero(n2)) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	mpz_tdiv_r_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
	return res;
}

integer& integer::operator %= (const integer& n)
{
	if (iszero(n)) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(*this)) return *this;
	integer res;   
	mpz_tdiv_r( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
	mpz_set( (mpz_ptr)&gmp_rep , (mpz_ptr)&(res.gmp_rep) );
	return *this;
}

integer& integer::operator %= (const unsigned long l)
{
	if (l ==0) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(*this)) return *this;
	mpz_tdiv_r_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
	return *this;
}

integer& integer::operator %= (const long l)
{
	if (l ==0) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(*this)) return *this;
	int sgn = SGN(l);
	mpz_tdiv_r_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, ABS(l));
	if (sgn <0) mpz_neg( (mpz_ptr)&gmp_rep, (mpz_ptr)&(gmp_rep) );
	return *this;
}


integer integer::operator % (const integer& n) const
{
	if (iszero(n)) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(*this)) return integer::zero;
	integer res;   
	mpz_tdiv_r( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
	return res;
}

long integer::operator % (const unsigned long l) const 
{
	if (l ==0) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(*this)) return 0L;
	integer Res(integer::one);   
	mpz_tdiv_r_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&gmp_rep, l);
	return integer2long( Res );
}

long integer::operator % (const long l) const 
{
	if (l ==0) {
		throw LinboxMathDivZero("[integer::/]: division by zero");
	}
	if (iszero(*this)) return 0L;
	integer Res(integer::one);   
	int sgn = SGN(l);
	mpz_tdiv_r_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&gmp_rep, ABS(l));
	if (sgn <0) Res = - Res;
	return integer2long( Res );
}
 
}
