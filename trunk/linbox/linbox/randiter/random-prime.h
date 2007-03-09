/* -*- mode:C++ -*- */
/*  File: primeiter.h
 *  Time-stamp: <09 Mar 07 16:29:32 Jean-Guillaume.Dumas@imag.fr> 
 */
#ifndef __RANDOM_PRIME_ITERATOR_H__
#define __RANDOM_PRIME_ITERATOR_H__
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
       
            /** @brief setSeed (unsigned long ul)
             *  Set the random seed to be ul.
             */
        void static setSeed(unsigned long ul) { 
	    integer::seeding(ul);
        }
        
        
    };
}

#endif
