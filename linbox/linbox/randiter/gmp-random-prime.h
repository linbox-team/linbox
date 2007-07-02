#ifndef __GMP_RANDOM_PRIME_H
#define __GMP_RANDOM_PRIME_H

#include <linbox/integer.h>

namespace LinBox {
/** \brief generating random prime integers, using the gmp library.

  * Similar to random-prime.h, but return type is integer, not long, allowing
  * for larger primes.
  * @author Dan Roche
  */

	class GmpRandomPrime {
	    public:
	    
		integer max;

		/** Create a new random prime generator.  Primes generated will
		  * be <= max.
		  */
		GmpRandomPrime( integer m ) :max(m) {}

		inline integer randomPrime() const {
			integer test;
			do test=test-1;//integer::nonzerorandom( test, max );
			while( !probab_prime( test, 10 ) );
			return test;
		}

		// I believe reference returned is appropriate. -bds
		//inline integer randomPrime( integer& p ) {
		inline integer& randomPrime( integer& p ) {
			do integer::nonzerorandom( p, max );
			while( !probab_prime( p, 10 ) );
			return p;
		}
		
	};
}

#endif
