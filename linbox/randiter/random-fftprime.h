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
#include "linbox/integer.h"
#include "linbox/util/error.h"
#include "linbox/util/timer.h"

namespace LinBox
{
	class RandomFFTPrime {
	public:
		typedef integer PrimeType; /* define the prime type */
		typedef std::vector<PrimeType> VectPrime;
		static const int32_t probab_prime_ntests = 25;
		static const size_t max_ntries = 256;

		/* Set p to an odd random prime such that
		 *  - p < pbound
		 *  - 2^k divides p-1
		 *
		 * The function repeatedly draws a random positive integer m and check
		 * if m*2^k+1 is prime, using Givaro::Protected::probab_prime.
		 *
		 * If after max_ntries tries, the function fails to find a prime, it
		 * returns false and p is set to 0.
		 */
		static inline bool randomPrime (PrimeType &p, const PrimeType &pbound, uint64_t k)
		{
			using Givaro::Protected::probab_prime;
			integer B = compute_max_m (pbound, k);
			if (B > 0)
			{
				for (size_t ntries = 0; ntries < max_ntries; ntries++)
				{
					integer::random_lessthan (p, B); /* 0 <= p < B */
					p = ((p+1) << k) + 1; /* 2^k+1 <= p < pbound */
					if (probab_prime (p, probab_prime_ntests))
						return true;
				}
			}
			p = 0;
			return false;
		}

		/* Set p to the largest odd prime satisfying
		 *  - p < pbound
		 *  - the 2-valuation of p-1 is maximal and >= kmin
		 *
		 * The function sets p to 0 and returns false, if there is no odd prime
		 * below pbound.
		 */
		static inline bool generatePrime (PrimeType &p, const PrimeType &pbound, uint64_t kmin = 1)
		{
			using Givaro::Protected::probab_prime;
			if (!kmin) /* set kmin to 1 if k == 0 */
				kmin++;
			for (uint64_t k = (pbound-2).bitsize()-1; k >= kmin; k--)
			{
				integer m = compute_max_m (pbound, k);
				/* test only odd value of m, so k is exactly the 2-valuation */
				if (!Givaro::isOdd (m))
					m--;
				for ( ; m > 0; m -= 2)
				{
					p = (m << k) + 1; /* p = m*2^k+1 */
					if (probab_prime (p, probab_prime_ntests))
						return true;
				}
			}
			p = 0;
			return false;
		}

		/* Set primes to a vector of distinct primes pi such that
		 *  - pi < pbound
		 *  - prod pi >= prodbound
		 *  - the 2-valuation of pi-1 is >= kmin and as large as possible
		 *
		 * The function returns false if there is no enough primes. In this
		 * case, the vector still contains the primes found so far.
		 */
		static inline bool generatePrimes (VectPrime &primes, const PrimeType &pbound, const PrimeType & prodbound, uint64_t kmin = 1)
		{
			using Givaro::Protected::probab_prime;
			primes.clear();
			PrimeType p, prod = 1;
			if (!kmin) /* set kmin to 1 if k == 0 */
				kmin++;
			for (uint64_t k = (pbound-2).bitsize()-1; k > kmin; k--)
			{
				integer m = compute_max_m (pbound, k);
				/* test only odd value of m, so k is exactly the 2-valuation */
				if (!Givaro::isOdd (m))
					m--;
				for ( ; m > 0; m -= 2)
				{
					p = (m << k) + 1; /* p = m*2^k+1 */
					if (probab_prime (p, probab_prime_ntests))
					{
						primes.push_back (p);
						prod *= p;
						if (prod > prodbound)
							return true;
					}
				}
			}
			return false;
		}

		/* Set the random seed */
		static inline void seeding (uint64_t seed)
		{
			integer::seeding (seed);
		}
		static inline void seeding (const Integer &seed)
		{
			integer::seeding (seed);
		}
		static inline void seeding ()
		{
			integer::seeding ();
		}

	private:
		RandomFFTPrime ()
		{
		}

		/* Given bound and k, compute maximum value of m such that
		 *      m*2^k+1 < bound  (Eq. 1)
		 *
		 * (Eq. 1) <=>   m < (bound-1)/2^k
		 *              <=>   m <= ((bound-1) >> k)-1      if 2^k | bound-1
		 *                    m <= (bound-1) >> k          otherwise
		 */
		static inline PrimeType compute_max_m (const PrimeType &bound, uint64_t k)
		{
			integer B = bound-1;
			B = (B & (uint64_t) ((1<<k) - 1)) ? (B >> k) : (B >> k) - 1;
			return B;
		}
	};
}

#endif //__LINBOX_random_fftprime_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
