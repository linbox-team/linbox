/* Copyright (C) 2010 LinBox
 *
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

#ifndef __LINBOX_gmp_random_prime_H
#define __LINBOX_gmp_random_prime_H

#include <linbox/integer.h>

namespace LinBox {
/** \brief generating random prime integers, using the gmp library.

  * Similar to random-prime.h, but return type is integer, not long, allowing
  * for larger primes.
  * @author Dan Roche
  */

	class GmpRandomPrime {
	    public:
	    
		integer max;

		/** Create a new random prime generator.  Primes generated will
		  * be <= max.
		  */
		GmpRandomPrime( integer m ) :max(m) {}

		inline integer randomPrime() const {
			integer test;
			do test=test-1;//integer::nonzerorandom( test, max );
			while( !probab_prime( test, 10 ) );
			return test;
		}

		// I believe reference returned is appropriate. -bds
		//inline integer randomPrime( integer& p ) {
		inline integer& randomPrime( integer& p ) {
			do integer::nonzerorandom( p, max );
			while( !probab_prime( p, 10 ) );
			return p;
		}
		
	};
}

#endif //__LINBOX_gmp_random_prime_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
