/** A class for generating random prime integers, using gmp library.
  * Similar to random-prime.h, but return type is integer, not long, allowing
  * for larger primes.
  * @author Dan Roche
  */

#ifndef __GMP_RANDOM_PRIME_H
#define __GMP_RANDOM_PRIME_H

#include <linbox/integer.h>

namespace LinBox {
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
			while( !mpz_probab_prime_p( test.get_rep(), 10 ) );
			return test;
		}

		inline integer randomPrime( integer& p ) {
			do integer::nonzerorandom( p, max );
			while( !mpz_probab_prime_p( p.get_rep(), 10 ) );
			return p;
		}
		
	};
}

#endif
