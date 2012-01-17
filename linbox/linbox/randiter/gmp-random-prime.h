/* Copyright (C) 2010 LinBox
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#error "deprecated and not tested"

#ifndef __LINBOX_gmp_random_prime_H
#define __LINBOX_gmp_random_prime_H

#include "linbox/integer.h"

namespace LinBox
{
	/** \brief generating random prime integers, using the gmp library.

	 * Similar to random-prime.h, but return type is integer, not long, allowing
	 * for larger primes.
	 * @author Dan Roche
	 * @deprecated This is not true.
	 */

	class GmpRandomPrime {
	public:

		integer max;

		/** Create a new random prime generator.  Primes generated will
		 * be <= max.
		 */
		GmpRandomPrime( integer m ) :
			max(m)
		{}

		inline integer randomPrime() const
		{
			integer test;
			do test=test-1;//integer::nonzerorandom( test, max );
			while( !probab_prime( test, 10 ) );
			return test;
		}

		// I believe reference returned is appropriate. -bds
		inline integer& randomPrime( integer& p )
		{
			do integer::nonzerorandom( p, max );
			while( !probab_prime( p, 10 ) );
			return p;
		}

	};
}

#endif //__LINBOX_gmp_random_prime_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

