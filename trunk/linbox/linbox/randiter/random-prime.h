/** -*- mode:C++ -*- */

/** File: RandomPrime.h
 *  Author: Zhendong Wan
 */

#ifndef __RANDOM_PRIME_H__
#define __RANDOM_PRIME_H__
#include <NTL/ZZ.h>

namespace LinBox {
	
	class RandomPrime {
		
	public:
		
		// define the prime type
		typedef long Prime_Type;

		/** @memo randomPrime()
		 *  return a random prime
		 */
		inline static Prime_Type randomPrime() {
			
			return NTL::GenPrime_long(30);
		}

		/** @memo randomPrime(Prime_Type& p)
		 *  return a random prime
		 */
		inline static Prime_Type randomPrime (Prime_Type& p) {

			return p = NTL::GenPrime_long(30);
			
		}
		
		/** @memo setSeed (unsigned long ul)
		 *  Set the random seed to be ul.
		 */
		void static setSeed(unsigned long ul) {
			
			NTL::SetSeed(NTL::to_ZZ(ul));
		}


	};
}

#endif
