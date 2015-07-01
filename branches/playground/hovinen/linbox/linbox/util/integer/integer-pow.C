/* -*- mode: c; style: linux -*- */

/* linbox/linbox/util/integer/integer-pow.C
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

Integer pow(const Integer& n, const unsigned long p)
{
	if (p == 0) return Integer::one;

//   Integer::Rep (res.gmp_rep)(l*ABS(n.gmp_rep.size));
	Integer Res;
	mpz_pow_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p);
	return Res;
#if 0
	int is_assg = 0;

	Integer res = Integer::one;
	Integer puiss  = n;

	while (p != 0) {
		if (p & 0x1) 
			if (is_assg) 
				res *= puiss;
			else { 
				is_assg = 1; 
				res = puiss; 
			}
		if ((p >>= 1) != 0) puiss = puiss * puiss;
	}
	return res;
#endif
}

Integer pow(const Integer& n, const long l) {
	if (l < 0)  return Integer::zero;
	return pow(n, (unsigned long) ABS(l) );
}


Integer powmod(const Integer& n, const unsigned long p, const Integer& m)
{
	if (p == 0) return Integer::one;
	Integer Res;
	mpz_powm_ui( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, p, (mpz_ptr)&m.gmp_rep);
	return Res;
}
Integer powmod(const Integer& n, const long e, const Integer& m)
{
	if (e < 0)  return Integer::zero;
	return powmod (n, (unsigned long)ABS(e), m);
}


Integer powmod(const Integer& n, const Integer& e, const Integer& m)
{
	if (e == 0) return Integer::one;
	if (e < 0)  return Integer::zero;
	Integer Res;
	mpz_powm( (mpz_ptr)&(Res.gmp_rep), (mpz_ptr)&n.gmp_rep, (mpz_ptr)&e.gmp_rep, (mpz_ptr)&m.gmp_rep);
	return Res;
}
 
}
