/* linbox/randiter/multimod-randomprime.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */


#ifndef __LINBOX_multimod_random_prime_H
#define __LINBOX_multimod_random_prime_H

#include "linbox/integer.h"
#include "linbox/util/timer.h"
#include <vector>

namespace LinBox
{

	class MultiModRandomPrime {
	protected:

		uint64_t   _bits;
		size_t   _size;
		integer _shift;

	public:


		MultiModRandomPrime(size_t n=1, uint64_t bits = 30, uint64_t seed = 0) :
			_bits(bits), _size(n)
		{
			_shift = integer(1)<<_bits;

			if (! seed)
				setSeed( (uint64_t) BaseTimer::seed() );
			else
				setSeed( seed );
		}

		// define the prime type
		typedef std::vector<integer> Prime_Type;

		/** \c randomPrime().
		 *  return a vector of random prime
		 */
		inline Prime_Type randomPrime() const
		{
			std::cout<<"ist used\n";
			std::vector<integer> tmp(_size);
			Givaro::IntPrimeDom IPD;
			for (size_t i=0;i<_size; ++i){
				do {
					integer::random(tmp[i],_bits-1);
					tmp[i]= _shift-tmp[i];
					IPD.nextprimein(tmp[i]);
			}
				while (std::find(tmp.begin(), tmp.begin()+(long)i, tmp[i]) != (tmp.begin()+(long)i ));
			}
			return tmp;
		}

		/** @c randomPrime(Prime_Type& p).
		 *  return a random prime
		 */
		inline Prime_Type randomPrime (Prime_Type& p) const
		{
			Givaro::IntPrimeDom IPD;
			for (size_t i=0;i<p.size(); ++i){
				do {
					integer::random(p[i],_bits-1);
					p[i]= _shift-p[i];
					IPD.nextprimein(p[i]);
			}
				while (std::find(p.begin(), p.begin()+(long)(i-1), p[i]) != (p.begin()+(long)(i-1) ));
			}
			return p;
		}


		/** @c setSeed(unsigned long ul).
		 *  Set the random seed to be ul.
		 */
		void static setSeed(unsigned long ul)
		{
			integer::seeding(ul);
		}


		inline Prime_Type createPrimes (const integer& maxPrime, const integer& productBound) const
		{
			Prime_Type tmp;
			integer acc=1, cur=maxPrime;
			Givaro::IntPrimeDom IPD;
			do {
				IPD.prevprimein(cur);
				tmp.push_back(cur);
				acc*=cur;
			} while (acc < productBound);

			return tmp;
		}



	};
}

#endif //__LINBOX_multimod_random_prime_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
