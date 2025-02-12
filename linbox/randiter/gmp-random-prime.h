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

#ifndef __LINBOX_gmp_random_prime_H
#define __LINBOX_gmp_random_prime_H

#include "linbox/randiter/random-prime.h"

namespace LinBox
{
	/** \brief generating random prime integers, using the gmp library.
	 * via primeiterator
	 */

	class GmpRandomPrime : public  PrimeIterator<> {
	public:
        GmpRandomPrime(const uint64_t bits = 23, const uint64_t seed = 0) :
                PrimeIterator(bits,seed) {}

		inline integer randomPrime()
		{
			this->operator++(); return this->operator*();
		}

		inline integer& randomPrime( integer& p )
		{
            return p=std::move(randomPrime());
		}

	};
}

#endif //__LINBOX_gmp_random_prime_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
