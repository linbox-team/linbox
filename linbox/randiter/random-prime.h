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
    namespace IteratorCategories {
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
    template<class IteratorTrait>
    struct UniqueSamplingTrait
		:public std::false_type { }; /* default: no guarantee of uniqueness */

    template<>
    struct UniqueSamplingTrait<IteratorCategories::DeterministicTag>
		:public std::true_type { };


        /*!  @brief  Prime Iterator.
         * @ingroup primes
         * @ingroup randiter
         *
         * Generates prime of specified length using a heuristically random distribution
         * (no guarantee whatsoever).
         * @internal
         * It is given by <code>nextprime(2^_bits-p)</code> where <code>size(p) < _bits</code>.
         */
    template<class Trait = IteratorCategories::HeuristicTag>
	class PrimeIterator{
	protected:
		uint64_t	_bits;  //!< common lenght of all primes
		integer        _prime;  //!< the generated prime.
		Givaro::IntPrimeDom _IPD; //!< empty struct dealing with primality.

        virtual void generatePrime();
	public:
		typedef integer Prime_Type ;
        typedef UniqueSamplingTrait<Trait> UniqueSamplingTag; //!< whether a prime can be picked more than once
        typedef Trait IteratorTag;

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
                _prime = Prime_Type(1)<<bits;
                generatePrime();
            }

            /** @brief operator++()  (prefix ++ operator)
             *  creates a new random prime.
             */
		inline PrimeIterator<Trait> &operator ++ ()
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

        uint64_t getBits() const {return _bits;}
	};

    template<>
    void PrimeIterator<IteratorCategories::HeuristicTag>::generatePrime(){
        integer::random_exact_2exp(_prime,_bits);
        _IPD.nextprimein(_prime);
        while (_prime.bitsize()>_bits)
            _IPD.prevprimein(_prime);
    }

    template<>
    void PrimeIterator<IteratorCategories::DeterministicTag>::generatePrime(){
        if (_prime < 3) {
            throw LinboxError("LinBox ERROR: Ran out of primes in Deterministic Prime Iterator.\n");
        }
        _IPD.prevprimein(_prime);
    }

    template<>
    void PrimeIterator<IteratorCategories::UniformTag>::generatePrime(){
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



/* ================
 * Find the log base 2 of an N-bit integer
 * From Bit Twiddling Hacks
 * By Sean Eron Anderson
 * seander@cs.stanford.edu
 * https://graphics.stanford.edu/~seander/bithacks.html
 */
    static const uint32_t MultiplyDeBruijnBitPosition[32] =
    {
        0U, 9U, 1U, 10U, 13U, 21U, 2U, 29U, 11U, 14U, 16U,
        18U, 22U, 25U, 3U, 30U, 8U, 12U, 20U, 28U, 15U, 17U,
        24U, 7U, 19U, 27U, 23U, 6U, 26U, 5U, 4U, 31U
    };

    uint32_t MultiplyDeBruijnHighestBit(uint32_t v)
    {
        v |= v >> 1; // first round down to one less than a power of 2
        v |= v >> 2;
        v |= v >> 4;
        v |= v >> 8;
        v |= v >> 16;

        return MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
    }
/* ================ */

	/*! @brief Adaptor class to make a fixed-length sequence behave like a PrimeIterator.
	 */
	template <typename IteratorT, class UniqueTrait = std::true_type>
	class PrimeSequence {
	private:
		IteratorT _cur;
		const IteratorT _end;

	public:
		using Prime_Type = typename std::remove_cv<typename std::remove_reference<decltype(*_cur)>::type>::type;
		using UniqueSamplingTag = UniqueTrait;
		using IteratorTag = IteratorCategories::DeterministicTag;

		PrimeSequence(IteratorT begin, IteratorT end) :
			_cur(begin), _end(end)
		{ }

		/** @brief operator++()  (prefix ++ operator)
		 *  moves to the next prime in the sequence.
		 */
		inline PrimeSequence& operator++ () {
			if (_cur == _end) {
				throw LinboxError("LinBox ERROR: Ran out of primes in PrimeSequence.\n");
			}
			++_cur;
			return *this;
		}

		/** @brief get the prime.
		 *  returns the actual prime.
		 *  @warning a new prime is not generated.
		 */
		const Prime_Type &operator * () const {return *_cur;}
	};

	/*! convenience factory to create a PrimeSequence from an STL-like container. */
	template <typename Container>
	PrimeSequence<typename Container::const_iterator>
		create_prime_sequence(const Container& cont)
	{
		return PrimeSequence<typename Container::const_iterator>(cont.begin(), cont.end());
	}


        /*!  @brief  Masked Prime Iterator.
         * @ingroup primes
         * @ingroup randiter
         *
         * Generates prime of specified length with fixed lower bits
         * @internal
         */
    template<class Trait = IteratorCategories::HeuristicTag>
    class MaskedPrimeIterator : public PrimeIterator<Trait> {
    private:
        const uint32_t  _shift;
        const uint32_t	_fffff;
        const uint32_t	_nmask;

    public:
        typedef MaskedPrimeIterator<Trait> Self_t;
        typedef PrimeIterator<Trait> Father_t;
        void generatePrime();

		inline Self_t &operator ++ ()
            {
                Father_t::operator++();
                return *this;
            }

            // Makes mask odd: implicit _mask=2*mask+1
        MaskedPrimeIterator(uint32_t mask, uint32_t max, uint64_t bits = 23, uint64_t seed = 0) :
                Father_t(bits,seed),
                _shift( MultiplyDeBruijnHighestBit(max) + 2),
                _fffff((1U<<_shift)-1),
                _nmask( ~((mask<<1) | 0x1) & _fffff ) // NOT _mask
            {
                this->_prime |= _fffff; // set lowest bits to 11111
                this->_prime ^= _nmask; // set lowest bits to _mask
                generatePrime();
            }

        const uint32_t getMask() const { return  ( ~(_nmask) & _fffff ); }
        const uint32_t getShift() const { return _shift; }
    };

    template<>
    void MaskedPrimeIterator<IteratorCategories::HeuristicTag>::generatePrime(){
        integer::random_exact_2exp(_prime,_bits);

        _prime |= _fffff; // set lowest bits to 11111
        _prime ^= _nmask; // set lowest bits to _mask

        while(! _IPD.isprime(_prime) ) {
            _prime += (1<<_shift);
        }
    }

    template<>
    void MaskedPrimeIterator<IteratorCategories::DeterministicTag>::generatePrime(){
        do {
            _prime -= (1<<_shift);
			if (_prime < 2) {
				throw LinboxError("LinBox ERROR: Ran out of primes in Masked Deterministic Prime Iterator.\n");
			}
        } while(! _IPD.isprime(_prime) );
    }


    template<>
    void MaskedPrimeIterator<IteratorCategories::UniformTag>::generatePrime(){
        do{
            integer::random_exact_2exp(_prime,_bits);
            switch (_prime %6){
                case 0: _prime++; break;
                case 4: _prime++; break;
                case 2: _prime--; break;
                case 3: _prime+=2; break;
            }

            _prime |= _fffff; // set lowest bits to 11111
            _prime ^= _nmask; // set lowest bits to _mask

        } while(!_IPD.isprime(_prime));
    }

   /*! @brief Adaptor class to make a single prime number behave like a PrimeIterator.
    */
   class FixedPrimeIterator {
    public:
        using Prime_Type = Givaro::Integer;
        using UniqueSamplingTag = std::false_type;
        using IteratorTag = IteratorCategories::DeterministicTag;

        FixedPrimeIterator (const Givaro::Integer& i) : _myprime(i) {}

        inline FixedPrimeIterator &operator ++ () { return *this; }

        const Prime_Type & operator * () const { return randomPrime(); }
        const Prime_Type & randomPrime() const { return _myprime; }

        void setBits(uint64_t bits) {}

        template<class _ModField> void setBitsField() {}

   private:
        const Prime_Type _myprime;
    };

}

#endif //__LINBOX_random_prime_iterator_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
