/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/util/integer/integer-cstor.C
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
#include "linbox/integer.h"

namespace LinBox
{

//------------------------------------- predefined null and one
const integer integer::zero(0UL);
const integer integer::one(1UL);


// -- integer(const char *s)
integer::integer(const char *s) 
{
	mpz_init_set_str((mpz_ptr)&gmp_rep, s, 10);
}


integer& integer::copy(const integer &n)
{
	if (this == &n) return *this;
	/* I don't understand what cpt is.  Suppose we just do the else clause?
	   if (cpt->decr() ==0) { 
	   mpz_clear((mpz_ptr)&gmp_rep) ; 
	   delete cpt ;
	   mpz_init_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
	   } else
	   -bds */
	mpz_set ( (mpz_ptr)&gmp_rep, (mpz_ptr)&(n.gmp_rep)) ;
	return *this ;
}

}
