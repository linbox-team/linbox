/* -*- mode:C++ -*- */

/*  File: RandomPrime.h
 *  Author: Zhendong Wan
 */

#ifndef __RANDOM_PRIME_H__
#define __RANDOM_PRIME_H__
#include <linbox/integer.h>
#include <linbox/util/timer.h>

namespace LinBox {
	
	class RandomPrime {		
	public:

	  int _bits;

            RandomPrime(int bits = 30, unsigned long seed = 0) 
                    : _bits(bits)
                {
                    if (! seed) 
                        RandomPrime::setSeed( BaseTimer::seed() );
                    else
                        RandomPrime::setSeed( seed );
                }
	  
	  // define the prime type
	  typedef integer Prime_Type;
	  
	  /** @memo randomPrime()
	   *  return a random prime
	   */
	  inline Prime_Type randomPrime() const {
	    integer tmp;
	    integer::random(tmp,_bits);
	    nextprime(tmp,tmp);
	    return tmp;
	  }
	  
	  /** @memo randomPrime(Prime_Type& p)
	   *  return a random prime
	   */
	  inline Prime_Type randomPrime (Prime_Type& p) const {
	    integer::random(p,_bits);
	    nextprime(p,p);
	    return p ; 
	  }
	  
	  /** @memo setSeed (unsigned long ul)
	   *  Set the random seed to be ul.
	   */
	  void static setSeed(unsigned long ul) { 
	    integer::seeding(ul);
	  }
	  
	  
	};
}

#endif
