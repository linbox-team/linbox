/*  -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen

// */
/*  File: primeiter.h
 *  Time-stamp: <12 Mar 07 18:35:26 Jean-Guillaume.Dumas@imag.fr> 
 */
#ifndef __LINBOX_random_prime_iterator_H
#define __LINBOX_random_prime_iterator_H
#include <linbox/integer.h>
#include <linbox/util/timer.h>

namespace LinBox {

	class RandomPrimeIterator {		
	public:

		int 	_bits;
		integer _shift;
		integer _prime;

		RandomPrimeIterator(int bits = 30, unsigned long seed = 0) 
		: _bits(bits)
		{
			_shift = integer(1)<<_bits;
			if (! seed) 
				RandomPrimeIterator::setSeed( BaseTimer::seed() );
			else
				RandomPrimeIterator::setSeed( seed );

			integer::random(_prime,_bits-1);
			_prime = _shift - _prime;
			nextprime( _prime, _prime);
		}

		RandomPrimeIterator(int low_bits = 30, int high_bits = 33,  unsigned long seed = 0) 
		{
			assert(low_bits<=high_bits);
			
			if (! seed) 
				RandomPrimeIterator::setSeed( BaseTimer::seed() );
			else
				RandomPrimeIterator::setSeed( seed );

			srand48(seed);
			double rnd = drand48();


			integer mp_low = 1 ;
			integer mp_hig = 1 ;

			mp_low << low_bits ;
			mp_hig << high_bits ;

			mp_diff = mp_hig-mp_low  ;



  { mpz_init(mp_rand); mpz_init(mp_temp); }
    { mpz_init(mp_lb);  mpz_init(mp_hb); }
      mpz_ui_pow_ui(mp_lb, 2, lb);
        mpz_ui_pow_ui(mp_hb, 2, hb);


			int rndbit

			_shift = integer(1)<<_bits;
			if (! seed) 
				RandomPrimeIterator::setSeed( BaseTimer::seed() );
			else
				RandomPrimeIterator::setSeed( seed );

			integer::random(_prime,_bits-1);
			_prime = _shift - _prime;
			nextprime( _prime, _prime);
		}

		// define the prime type
		typedef integer Prime_Type;

		/** @brief operator++()
		 *  creates a new random prime
		 */
		inline RandomPrimeIterator &operator ++ () { 
			integer::random(_prime,_bits-1);
			_prime = _shift - _prime;
			nextprime( _prime, _prime);
			return *this;
		}

		/** @brief operator*()
		 *  returns the actual prime
		 */
		Prime_Type &operator *  () { 
			return _prime;
		}
		Prime_Type & randomPrime() { 
			return _prime;
		}

		/** @brief setSeed (unsigned long ul)
		 *  Set the random seed to be ul.
		 */
		void static setSeed(unsigned long ul) { 
			integer::seeding(ul);
		}


	};
}

#endif //__LINBOX_random_prime_iterator_H
