/* Copyright (C) 2007,2010 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr>
 * Modified by Brice Boyer (briceboyer) <boyer.brice@gmail.com> (RandomPrimeIter)
 * Modified by Clement Pernet <clement.pernet@imag.fr> (PrimeIter)
 *
 *
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

/*! @file randiter/random-prime.h
 * @ingroup randiter
 * @brief Generates random positive prime \ref integers.
 */

#ifndef __LINBOX_random_prime_iterator_H
#define __LINBOX_random_prime_iterator_H


#include <givaro/givinteger.h>
#include <givaro/givintprime.h>
#include "linbox/util/timer.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"
namespace LinBox
{
            /*! \brief Information about the type of Prime Iterator
             */
        namespace RandomCategories {
                //! Iterator following a deterministic sequence of primes (from the largest one, in decreasing order
                struct DeterministicTag{};
                //! Iterator sampling uniformly from all primes of given bitsize
                struct UniformTag{};
                //! Iterator sampling randomly (no distribution guaranteed whatsoever) from all primes of given bitsize
                struct HeuristicTag{};
        };
        
            /*! \brief Whether a prime generator generates a sequence with non repeating
             * numbers
             */
        template<class RandomTrait>
        struct UniqueSamplingTrait;

        template<>
        struct UniqueSamplingTrait<RandomCategories::DeterministicTag>{
                typedef std::true_type value;
        };
        template<>
        struct UniqueSamplingTrait<RandomCategories::UniformTag>{
                typedef std::false_type value;
        };
        template<>
        struct UniqueSamplingTrait<RandomCategories::HeuristicTag>{
                typedef std::false_type value;
        };

        
        /*!  @brief  Prime Iterator.
	 * @ingroup primes
	 * @ingroup randiter
	 *
	 * Generates prime of specified length using a heuristically random distribution 
         * (no guarantee whatsoever).
	 * @internal
	 * It is given by <code>nextprime(2^_bits-p)</code> where <code>size(p) < _bits</code>.
	 */
        template<class RandomTrait = RandomCategories::HeuristicTag>
	class PrimeIterator{
	private:
		uint64_t 	_bits;  //!< common lenght of all primes
		integer        _prime;  //!< the generated prime.
		Givaro::IntPrimeDom _IPD; //!< empty struct dealing with primality.

                void generatePrime();
	public:
		typedef integer Prime_Type ;
                typedef typename UniqueSamplingTrait<RandomTrait>::value UniqueSamplingTag; //!< whether a prime can be picked more than once
                typedef RandomTrait RandomTrait_Type;

                /*! Constructor.
		 * @param bits size of primes (in bits). Default is 23 so it
		 * can fit in a <code>Linbox::Modular<double></code>.
		 * @param seed if \c 0 a seed will be generated, otherwise, the
		 * provided seed will be use.
		 */
		PrimeIterator(uint64_t bits = 23, uint64_t seed = 0) :
			_bits(bits)
		{
			linbox_check(bits >1);
			if (! seed)
                                seed = BaseTimer::seed();
                        setSeed (seed);
                        generatePrime();
		}

                /** @brief operator++()  (prefix ++ operator)
		 *  creates a new random prime.
		 */
		inline PrimeIterator<RandomCategories::HeuristicTag> &operator ++ ()
		{
			generatePrime();
                        return *this;
		}

		/** @brief get the random prime.
		 *  returns the actual prime.
		 *  @warning a new prime is not generated.
		 */
		const Prime_Type &operator * () const {return _prime;}

                /** @brief Sets the seed.
		 *  Set the random seed to be \p ul.
		 *  @param ul the new seed.
		 */
		void static setSeed(uint64_t ul){integer::seeding(ul);}

                /** @brief Sets the bit size.
		 *  @param bits the new bit size.
		 */
		void setBits(uint64_t bits) {
			_bits = bits;
			linbox_check(bits >1);
			generatePrime();
		}

                /** @brief Sets the bit size.
		 *  @param bits the new bit size.
		 */
		template<class _ModField>
		void setBitsField()
		{
			integer k = FieldTraits<_ModField >::maxModulus();
			uint64_t bits = (uint64_t)(k.bitsize());
			if (!bits)
				throw("weird");
			--bits;
			if (bits < _bits)
				setBits(bits);
		}

		template<class _ModField>
		bool setBitsDelayedField(size_t n)
		{
			integer k = FieldTraits<_ModField >::maxModulus();
			int bits = (int)(k.bitsize());
			if (!bits) throw("weird");
			--bits;
			bits -= (int)((log((double)n)/2./M_LN2));
			if (bits < 0) return false;
			if (bits < (int64_t)_bits) setBits((uint64_t)bits);
			return true;
		}
	};

        template<>
        void PrimeIterator<RandomCategories::HeuristicTag>::generatePrime(){
                integer::random_exact_2exp(_prime,_bits);
                _IPD.nextprimein(_prime);
                while (_prime.bitsize()>_bits)
                        _IPD.prevprimein(_prime);
        }

        template<>
        void PrimeIterator<RandomCategories::DeterministicTag>::generatePrime(){_IPD.prevprimein(_prime);}

        template<>
        void PrimeIterator<RandomCategories::UniformTag>::generatePrime(){
                do{
                        integer::random_exact_2exp(_prime,_bits);
                        switch (_prime %6){
                            case 0: _prime++; break;
                            case 4: _prime++; break;
                            case 2: _prime--; break;
                            case 3: _prime+=2; break;
                        }
                } while(!_IPD.isprime(_prime));
        }
}

#endif //__LINBOX_random_prime_iterator_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
