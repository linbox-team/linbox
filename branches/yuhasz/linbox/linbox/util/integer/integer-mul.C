/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/util/integer/integer-mul.C
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

namespace LinBox
{

//-------------------------------------------------- operator *
integer& integer::mulin(integer& res, const integer& n) 
{
	if (iszero(n)) return res = integer::zero;
	if (iszero(res)) return res;
	mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
	return res;
}
integer& integer::mulin(integer& res, const long n) 
{
	if (iszero(n)) return res = integer::zero;
	if (iszero(res)) return res;
//   int sgn = SGN(n); 
	int sgn = SGN(n); 
	mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, ABS(n));
//   if (sgn <0) res.gmp_rep.size = -res.gmp_rep.size;
	if (sgn <0) return res = -res;
	return res;
}
integer& integer::mulin(integer& res, const unsigned long n) 
{
	if (iszero(n)) return res = integer::zero;
	if (iszero(res)) return res;
	mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
	return res;
}

integer& integer::mul(integer& res, const integer& n1, const integer& n2)
{
	if (iszero(n1)) return res = integer::zero;
	if (iszero(n2)) return res = integer::zero;
	mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
	return res;
}
integer& integer::mul(integer& res, const integer& n1, const long n2)
{
	if (iszero(n1)) return res = integer::zero;
	if (iszero(n2)) return res = integer::zero;
	int sgn = SGN(n2); 
	mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, ABS(n2));
//   if (sgn <0) res.gmp_rep.size = -res.gmp_rep.size;
	if (sgn <0) return res = -res;
	return res;
}
integer& integer::mul(integer& res, const integer& n1, const unsigned long n2)
{
	if (iszero(n1)) return res = integer::zero;
	if (iszero(n2)) return res = integer::zero;
	mpz_mul_ui( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
	return res;
}

integer& integer::axpy(integer& res, const integer& a, const integer& x, const integer& b)
{
	if (&res == &b) return integer::axpyin(res,a,x);
	if (iszero(a) || iszero(x)) return res = b;
	mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&a.gmp_rep, (mpz_ptr)&x.gmp_rep);
	mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&b.gmp_rep);
	return res;
}

integer& integer::axpyin(integer& res, const integer& a, const integer& x)
{
	if (iszero(a) || iszero(x)) return res;
	Rep gmp_res; mpz_init((mpz_ptr)&gmp_res);
	mpz_mul( (mpz_ptr)&gmp_res, (mpz_ptr)&a.gmp_rep, (mpz_ptr)&x.gmp_rep);
	mpz_add( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&gmp_res);
	mpz_clear((mpz_ptr)&gmp_res);
	return res;
}


integer& integer::axmy(integer& res, const integer& a, const integer& x, const integer& b)
{
	if (&res == &b) return integer::axmyin(res,a,x);
	if (iszero(a) || iszero(x)) return integer::neg(res,b);
	mpz_mul( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&a.gmp_rep, (mpz_ptr)&x.gmp_rep);
	mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&b.gmp_rep);
	return res;
}
integer& integer::axmyin(integer& res, const integer& a, const integer& x)
{
	if (iszero(a) || iszero(x)) return res;
	Rep gmp_res; mpz_init((mpz_ptr)&gmp_res);
	mpz_mul( (mpz_ptr)&gmp_res, (mpz_ptr)&a.gmp_rep, (mpz_ptr)&x.gmp_rep);
	mpz_sub( (mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&gmp_res);
	mpz_clear((mpz_ptr)&gmp_res);
	return res;
}



integer& integer::operator *= (const integer& n)
{
	if (iszero(n)) return *this = integer::zero;
	if (iszero(*this)) return *this;
//   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );   
	integer res; 
	mpz_mul( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
	return *this = res;
}

integer& integer::operator *= (const unsigned long l)
{
	if (l==0) return *this = integer::zero;
	if (iszero(*this)) return *this;
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
	mpz_mul_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, l);
	return *this;
}

integer& integer::operator *= (const long l)
{
	if (l==0) return *this =integer::zero;
	if (iszero(*this)) return *this;
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );
	int sgn = SGN(l);
	mpz_mul_ui( (mpz_ptr)&(gmp_rep), (mpz_ptr)&gmp_rep, ABS(l));
	if (sgn <0) mpz_neg( (mpz_ptr)&gmp_rep, (mpz_ptr)&(gmp_rep) );
	return *this;
}


integer integer::operator * (const integer& n) const
{
	if (iszero(n)) return integer::zero;
	if (iszero(*this)) return integer::zero;
//   Rep (res.gmp_rep)( MAX(SZ_REP(n.gmp_rep),SZ_REP(gmp_rep)) );   
	integer res;
	mpz_mul( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
	return res;
}

integer integer::operator * (const unsigned long l) const 
{
	if (l==0) return integer::zero;
	if (iszero(*this)) return integer::zero;
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
	integer res;
	mpz_mul_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, l);
	return res;
}

integer integer::operator * (const long l) const 
{
	if (l==0) return integer::zero;
	if (iszero(*this)) return integer::zero;
//   Rep (res.gmp_rep)( MAX(SZ_REP(gmp_rep),1) );   
	integer res;   
	int sgn = SGN(l);
	mpz_mul_ui( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&gmp_rep, ABS(l));
//   if (sgn <0) (res.gmp_rep).size = -(res.gmp_rep).size;
//   return integer((res.gmp_rep));
	if (sgn <0) mpz_neg( (mpz_ptr)&(res.gmp_rep), (mpz_ptr)&(res.gmp_rep) );
	return res;
}
 
}
