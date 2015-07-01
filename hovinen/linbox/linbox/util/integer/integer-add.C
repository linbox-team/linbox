/* -*- mode: c; style: linux -*- */

/* linbox/linbox/util/integer/integer-add.C
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

#include <gmp.h>

#include "linbox/integer.h"

namespace LinBox 
{

Integer& Integer::addin (Integer& res, const Integer& n) 
{
	if (iszero (n)) return res;
	if (iszero (res)) return res = n;
	mpz_add ((mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, (mpz_ptr)&n.gmp_rep );
	return res;
}

Integer& Integer::addin (Integer& res, const long n) 
{
	if (iszero (n)) return res;
	if (iszero (res)) return res = n;
	int sgn = SGN (n); 
	if (sgn >0) mpz_add_ui ((mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, n);
	else mpz_sub_ui ((mpz_ptr)&res.gmp_rep, (mpz_ptr)&res.gmp_rep, -n);
	return res;
}

Integer& Integer::addin (Integer& res, const unsigned long n) 
{
	if (iszero (n)) return res;
	if (iszero (res)) return res = n;
	mpz_add_ui ((mpz_ptr)&res.gmp_rep, (mpz_srcptr)&res.gmp_rep, n);
	return res;
}

Integer& Integer::add (Integer& res, const Integer& n1, const Integer& n2)
{
	if (iszero (n1)) return res = n2;
	if (iszero (n2)) return res = n1;
	mpz_add ((mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, (mpz_ptr)&n2.gmp_rep);
	return res;
}

Integer& Integer::add (Integer& res, const Integer& n1, const long n2)
{
	if (iszero (n1)) return res = n2;
	if (iszero (n2)) return res = n1;
	int sgn = SGN (n2); 
	if (sgn >0) mpz_add_ui ((mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
	else mpz_sub_ui ((mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, -n2);
	return res;
}

Integer& Integer::add (Integer& res, const Integer& n1, const unsigned long n2)
{
	if (iszero (n1)) return res = n2;
	if (iszero (n2)) return res = n1;
	mpz_add_ui ((mpz_ptr)&res.gmp_rep, (mpz_ptr)&n1.gmp_rep, n2);
	return res;
}

Integer& Integer::operator += (const Integer& n)
{
	if (iszero (n)) return *this;
	if (iszero (*this)) return logcpy (n);
	mpz_add ((mpz_ptr) &(gmp_rep), (mpz_ptr)&gmp_rep, (mpz_ptr)&n.gmp_rep) ;
	return *this;
}

Integer& Integer::operator += (const unsigned long l)
{
	if (l==0) return *this;
	if (iszero (*this)) return logcpy (Integer (l));
	mpz_add_ui ((mpz_ptr) &(gmp_rep), (mpz_ptr)&gmp_rep, l);
	return *this;
}

Integer& Integer::operator += (const long l)
{
	if (l==0) return *this;
	if (iszero (*this)) return logcpy (Integer (l));
	int sgn = SGN (l);
	if (sgn >0) mpz_add_ui ((mpz_ptr) &(gmp_rep), (mpz_ptr)&gmp_rep, l);
	else mpz_sub_ui ((mpz_ptr) &(gmp_rep), (mpz_ptr)&gmp_rep, -l);
	return *this;
}


Integer Integer::operator + (const Integer& n) const
{
	if (iszero (n)) return *this;
	if (iszero (*this)) return n;
	Integer res;
	mpz_add ((mpz_ptr) &(res.gmp_rep), (mpz_ptr) &gmp_rep, (mpz_ptr) &n.gmp_rep);
	return res;
}

Integer Integer::operator + (const unsigned long l) const 
{
	if (l == 0) return *this;
	if (iszero (*this)) return Integer (l);
	Integer res;   
	mpz_add_ui ((mpz_ptr) &(res.gmp_rep), (mpz_ptr) &gmp_rep, l);
	return res;
}

Integer Integer::operator + (const long l) const 
{
	if (l == 0) return *this;
	if (iszero (*this)) return Integer (l);
	Integer res;   
	int sgn = SGN (l);
	if (sgn >0) mpz_add_ui ((mpz_ptr) &(res.gmp_rep), (mpz_ptr) &gmp_rep, l);
	else mpz_sub_ui ((mpz_ptr) &(res.gmp_rep), (mpz_ptr) &gmp_rep, -l);
	return res;
}
 
}
