/* linbox/randiter/multimod-randomprime.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */


#ifndef __LINBOX_multimod_random_prime_H
#define __LINBOX_multimod_random_prime_H

#include <linbox/integer.h>
#include <linbox/util/timer.h>
#include <vector>

namespace LinBox {
	
	class MultiModRandomPrime {		
	protected:

		size_t   _bits;
		size_t   _size;
		integer _shift;

	public:
		

		MultiModRandomPrime(size_t n=1, size_t bits = 30, unsigned long seed = 0) 
			: _bits(bits), _size(n)
		{
			_shift = integer(1)<<_bits;

			if (! seed) 
				MultiModRandomPrime::setSeed( BaseTimer::seed() );
			else
				MultiModRandomPrime::setSeed( seed );
		}
	  
		// define the prime type
		typedef std::vector<integer> Prime_Type;
	  
		/** @memo randomPrime()
		 *  return a vector of random prime
		 */
		inline Prime_Type randomPrime() const {std::cout<<"ist used\n";
			std::vector<integer> tmp(_size);
			for (size_t i=0;i<_size; ++i){
				do {
					integer::random(tmp[i],_bits-1);
					tmp[i]= _shift-tmp[i];
					nextprime(tmp[i],tmp[i]);					
				} while (find(tmp.begin(), tmp.begin()+i, tmp[i]) != (tmp.begin()+i ));				
			}				
			return tmp;
		}
	  
		/** @memo randomPrime(Prime_Type& p)
		 *  return a random prime
		 */
		inline Prime_Type randomPrime (Prime_Type& p) const {
			for (size_t i=0;i<p.size(); ++i){
				do {
					integer::random(p[i],_bits-1);
					p[i]= _shift-p[i];
					nextprime(p[i],p[i]);
				} while (find(p.begin(), p.begin()+i-1, p[i]) != (p.begin()+i-1 ));		
			}
			return p;
		}
	  

		/** @memo setSeed (unsigned long ul)
		 *  Set the random seed to be ul.
		 */
		void static setSeed(unsigned long ul) { 
			integer::seeding(ul);
		}


		inline Prime_Type createPrimes (const integer& maxPrime, const integer& productBound) const {
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

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
