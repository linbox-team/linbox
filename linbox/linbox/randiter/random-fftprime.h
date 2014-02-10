/* linbox/algorithms/
 * Copyright (C) 2005  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
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


#ifndef __LINBOX_random_fftprime_H
#define __LINBOX_random_fftprime_H
#include <vector>
using namespace std;
#include "linbox/integer.h"
#include "linbox/util/debug.h"
#include "linbox/util/timer.h"

namespace LinBox
{

	class RandomFFTPrime {
	public:

		size_t         _bits;

		RandomFFTPrime(size_t bits = 20, unsigned long seed = 0) :
			_bits(bits)
		{
			if (! seed)
				RandomFFTPrime::setSeed( (unsigned long)BaseTimer::seed() );
			else
				RandomFFTPrime::setSeed( seed );
		}

		// define the prime type
		typedef integer Prime_Type;

                 /** @brief randomPrime(size_t b)
		 *  return a random FFT prime with a 2-valuation larger than b in its order
                 *  the randomness is on the FFT primes lying in the given range
                 *  an error is thrown if no such prime exist
		 */
		inline Prime_Type randomPrime (size_t b) const
		{
                        integer tmp;
                        randomPrime(tmp,b);
			return tmp;
                }

                /** @brief randomPrime(Prime_Type& p, size_t b)
		 *  return a random FFT prime with a 2-valuation larger than b in its order
                 *  the randomness is on the FFT primes lying in the given range
                 *  an error is thrown if no such prime exist
		 */
		inline Prime_Type randomPrime (Prime_Type& t, size_t b) const
		{
                        linbox_check(b<_bits);
			size_t tresh;
			do {
                                size_t cbits= (size_t)rand() %(_bits-b);
				tresh = 1<<(cbits);
				size_t p = 1<<((size_t)_bits-cbits);
				do {
					integer::random(t,cbits);
					t = t*p+1;
					tresh--;
				} while (!Givaro::probab_prime(t,25) && (tresh));
			}
			while(tresh==0);
                        linbox_check(Givaro::probab_prime(t,25))
			return t;
		}

		/** @brief generatePrime()
		 *  return a FFT prime with the largest 2-valuation in its order
		 */
		inline Prime_Type generatePrime() const
		{
			integer tmp;
                        generatePrime(tmp);
			return tmp;
		}

		/** @brief generatePrime(Prime_Type& p)
		 *  return a FFT prime with the largest 2-valuation in its order
		 */
		inline Prime_Type generatePrime (Prime_Type& t) const
		{
			size_t cbits=5;
			size_t tresh;
			do {
				tresh = 1<<(cbits);
				size_t p = 1<<((size_t)_bits-cbits);
				do {
					integer::random(t,cbits);
					t = t*p+1;
					tresh--;
				} while (!Givaro::probab_prime(t,25) && (tresh));
				cbits++;
			}
			while(tresh==0);

			return t;
		}

                // generate a vector of distinct FFT primes with largest 2-valuation
                inline std::vector<Prime_Type> generatePrimes (std::vector<Prime_Type>& primes) const {
                        size_t pos = 0;
                        size_t k= primes.size();
                        integer tmp;
                        for (long b = (long)_bits - 1; b >= 0; b--)
                                for (long l = (1L << ((long)_bits - b - 1)) + 1; l < (1L << ((long)_bits - b)); l +=2) {
                                        tmp = (1L << b) * l + 1;
                                        if (Givaro::probab_prime(tmp, 25) >= 1) {
                                                primes[pos] = tmp;
                                                pos++;
                                                if (pos >= k)
                                                        return primes;
                                        }
                                }
                        linbox_check(primes[k] != 0); // Could not find enough primes
                        return primes;
                }

		/** @brief setSeed (unsigned long ul)
		 *  Set the random seed to be ul.
		 */
		void static setSeed(unsigned long ul)
		{
			integer::seeding(ul);
		}
	};
}

#endif //__LINBOX_random_fftprime_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

