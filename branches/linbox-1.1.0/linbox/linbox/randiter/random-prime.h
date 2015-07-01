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
	  integer _shift;

            RandomPrime(int bits = 30, unsigned long seed = 0) 
                    : _bits(bits)
                {
		  _shift = integer(1)<<_bits;
                    if (! seed) 
                        RandomPrime::setSeed( BaseTimer::seed() );
                    else
                        RandomPrime::setSeed( seed );
                }
	  
	  // define the prime type
	  typedef integer Prime_Type;
	  
	  /** @brief randomPrime()
	   *  return a random prime
	   */
	  inline Prime_Type randomPrime() const {
	    integer tmp;
	    integer::random(tmp,_bits-1);
	    tmp= _shift - tmp;
	    nextprime(tmp,tmp);
	    return tmp;
	  }
	  
	  /** @brief randomPrime(Prime_Type& p)
	   *  return a random prime
	   */
	  inline Prime_Type randomPrime (Prime_Type& p) const {
	    integer::random(p,_bits-1);
	    p= _shift - p;
	    nextprime(p,p);
	    return p ; 
	  }
	  
	  /** @brief setSeed (unsigned long ul)
	   *  Set the random seed to be ul.
	   */
	  void static setSeed(unsigned long ul) { 
	    integer::seeding(ul);
	  }

	  	  
	};
}

#endif
