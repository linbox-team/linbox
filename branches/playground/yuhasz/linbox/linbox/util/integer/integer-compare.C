/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/util/integer/integer-compare.C
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

int compare (const integer &a, const integer& b) 
{
	if (&a == &b) return 0 ;

	int cmp = mpz_cmp ( (mpz_ptr)&a.gmp_rep, (mpz_ptr)&b.gmp_rep );
	return (cmp <0 ? -1 : (cmp >0 ? 1 : 0));
}

int absCompare (const integer &a, const integer &b) 
{
	integer c = ((a<0)?-a:a), d = ((b<0)?-b:b);
	return mpz_cmp ( (mpz_ptr)&(c.gmp_rep), (mpz_ptr)&(d.gmp_rep));
}

int integer::operator != (const int l) const 
{ return mpz_cmp_si ( (mpz_ptr)&gmp_rep, l ) != 0; }

int integer::operator != (const long l) const 
{ return mpz_cmp_si ( (mpz_ptr)&gmp_rep, l ) != 0; }

int integer::operator > (const int l) const 
{ return mpz_cmp_si ((mpz_ptr)&gmp_rep, l) > 0; }

int integer::operator > (const long l) const 
{ return mpz_cmp_si ((mpz_ptr)&gmp_rep, l) > 0; }

int integer::operator < (const int l) const 
{ return mpz_cmp_si ((mpz_ptr)&gmp_rep, l) < 0; }

int integer::operator < (const long l) const 
{ return mpz_cmp_si ((mpz_ptr)&gmp_rep, l) < 0; }
 
}
