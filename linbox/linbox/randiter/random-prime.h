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

	  int _bits;

	  RandomPrime(int bits = 30)
	  {
	    _bits = bits;
	  };
	  
	  // define the prime type
	  typedef long Prime_Type;
	  
	  /** @memo randomPrime()
	   *  return a random prime
	   */
	  inline Prime_Type randomPrime() const {
	    return NTL::GenPrime_long(_bits);
	  }
	  
	  /** @memo randomPrime(Prime_Type& p)
	   *  return a random prime
	   */
	  inline Prime_Type randomPrime (Prime_Type& p) const {
	    
	    return p = NTL::GenPrime_long(_bits);
	    
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
