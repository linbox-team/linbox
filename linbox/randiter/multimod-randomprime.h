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

		size_t   _bits;
		size_t   _size;
		integer _shift;

	public:


		MultiModRandomPrime(size_t n=1, size_t bits = 30, unsigned long seed = 0) :
			_bits(bits), _size(n)
		{
			_shift = integer(1)<<_bits;

			if (! seed)
				MultiModRandomPrime::setSeed( (unsigned long) BaseTimer::seed() );
			else
				MultiModRandomPrime::setSeed( seed );
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
			for (size_t i=0;i<_size; ++i){
				do {
					integer::random(tmp[i],_bits-1);
					tmp[i]= _shift-tmp[i];
					nextprime(tmp[i],tmp[i]);
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
			for (size_t i=0;i<p.size(); ++i){
				do {
					integer::random(p[i],_bits-1);
					p[i]= _shift-p[i];
					nextprime(p[i],p[i]);
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

			do {
				prevprime(cur, cur);
				tmp.push_back(cur);
				acc*=cur;
			} while (acc < productBound);

			return tmp;
		}



	};
}

#endif //__LINBOX_multimod_random_prime_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

