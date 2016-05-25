/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
#include "linbox/util/debug.h"
#include "linbox/util/timer.h"

namespace LinBox
{

	class RandomFFTPrime {
	public:
		// define the prime type
		typedef integer Prime_Type;

		uint64_t           _bits;
		Prime_Type  _prime_bound;

		RandomFFTPrime(Prime_Type pbound=0x100000, unsigned long seed = 0) :
			_bits(pbound.bitsize()), _prime_bound(pbound)
		{
			if (! seed)
				RandomFFTPrime::setSeed( (unsigned long)BaseTimer::seed() );
			else
				RandomFFTPrime::setSeed( seed );
		}

		
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
		inline Prime_Type randomPrime (Prime_Type& t, uint64_t b) const
		{
			linbox_check(b<_bits);
			size_t tresh;
			do {
				size_t cbits= (size_t)rand() %(_bits-b);
				tresh = 1<<(cbits);
				uint64_t p = 1<<((size_t)_bits-cbits);
				do {
					integer::random(t,cbits);
					t = t*integer(p)+1;
					tresh--;
				} while (!Givaro::Protected::probab_prime(t,25) && (tresh));
			}
			while(tresh==0);
			linbox_check(Givaro::Protected::probab_prime(t,25))
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
			size_t cbits=1;
			size_t tresh;
			do {
				tresh = 1<<(cbits);
				uint64_t p = 1<<((size_t)_bits-cbits);
				do {
					integer::random(t,cbits);
					t = t*integer(p)+1;
					tresh--;
				} while (!Givaro::Protected::probab_prime(t,25) && (tresh));
				cbits++;
			}
			while(tresh==0);

			return t;
		}

		// generate a vector of distinct FFT primes with largest 2-valuation
		// s.t. their product is larger than a given bound
		inline std::vector<Prime_Type> generatePrimes (const Prime_Type & bound) const {
			std::vector<Prime_Type> primes;
			Prime_Type prod=1;
			integer tmp;
			for (int64_t b = _bits - 1; b >= 0; b--)
				for (int64_t l = ((int64_t)1 << (_bits - b - 1)) + 1; l < (1L << (_bits - b)); l +=2) {
					tmp = ((int64_t)1 << b) * l + 1;
					if (Givaro::Protected::probab_prime(tmp, 25) >= 1) {
						primes.push_back(tmp);
						prod*=tmp;
						if (prod > bound)
							return primes;
					}
				}
			linbox_check(prod > bound ); // Could not find enough primes
			return primes;
		}

		// generate a vector of distinct FFT primes with largest 2-valuation
		// s.t. their product is larger than a given bound
		inline bool generatePrimes (const Prime_Type & bound, std::vector<Prime_Type> &primes) const {
			primes.clear();
			Prime_Type prod=1;
			integer tmp;
			for (int64_t b = (int64_t)_bits - 1; b >= 0; b--)
				for (int64_t l = (1L << ((int64_t)_bits - b - 1)) + 1; l < (1L << ((int64_t)_bits - b)); l +=2) {
					tmp = (1L << b) * l + 1;
					if (Givaro::Protected::probab_prime(tmp, 25) >= 1) {
						primes.push_back(tmp);
						prod*=tmp;
						if (prod > bound){
							return true;
						}
					}
				}
			return false; // false -> Could not find enough primes
		}

		size_t twoVal(integer t) const {
			integer x=t;
			size_t v=0;
			while(x%2 == 0) {v++;x/=2;}
			return v;
		}

		// generate a vector of distinct FFT primes with  2-valuation largest than val
		// s.t. their product is larger than a given bound
		inline bool generatePrimes ( uint64_t val, const Prime_Type & bound, std::vector<Prime_Type> &primes) const {
			primes.clear();
			Prime_Type prod=1;
			integer tmp;
			// std::cout<<"rns bound: "<<bound<<std::endl;
			// std::cout<<"2 valuation: "<<val<<std::endl;
			// std::cout<<"prime bitmax: "<<_bits<<std::endl;
			// std::cout<<"prime max: "<<_prime_bound<<std::endl;

			if (val > _bits) return false;

#if 0
			for (int64_t b = (int64_t)_bits; b >= (int64_t)val; b--)
				// for (uint64_t l = (1ULL << ((int64_t)_bits - b - 1)) + 1; l < (1ULL << ((int64_t)_bits - b)); l +=2) {
				for (int64_t l = ((int64_t)1 << ((int64_t)_bits - b)) - 1; l >=1; l -=2) {
					tmp = ((int64_t)1 << b) * l + 1;
					if (Givaro::Protected::probab_prime(tmp, 25) >= 1) {
						primes.push_back(tmp);
						prod*=tmp;
						//std::cout<<tmp<<" -> "<<tmp.bitsize()<<" (order="<<twoVal(tmp-1)<<") "<<prod<<std::endl;
						if (prod > bound){
							return true;
						}
					}
				}
#else
			for (int64_t l = (_prime_bound -1) >>val ; l >=1; l -=1) {
				tmp = ((int64_t)1 << val) * l + 1;
				if (Givaro::Protected::probab_prime(tmp, 25) >= 1) {
					primes.push_back(tmp);
					prod*=tmp;
					//std::cout<<tmp<<" -> "<<tmp.bitsize()<<" (order="<<twoVal(tmp-1)<<") "<<prod<<std::endl;
					if (prod > bound){
						// try to replace the last prime with a smallest one
						for (int64_t k=1;k<l;k++){
							tmp = ((int64_t)1 << val) * k + 1;
							if (Givaro::Protected::probab_prime(tmp, 25) >= 1) {
								if (prod*tmp > bound*primes.back()){
									//std::cout<<"replacing prime "<<primes.back()<<" with "<<tmp<< " -> "<<tmp.bitsize()<<" (order="<<twoVal(tmp-1)<<") ";
									prod/=primes.back();
									primes.back()=tmp;
									prod*=tmp;
									//std::cout<<prod<<std::endl;
									return true;
								}
							}
						}

						return true;
					}
				}
			}


#endif
			return false; // false -> Could not find enough primes
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
