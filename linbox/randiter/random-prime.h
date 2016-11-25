/* Copyright (C) 2007,2010 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr>
 * Modified by Brice Boyer (briceboyer) <boyer.brice@gmail.com> (RandomPrimeIter)
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

	/*!  @brief Random Prime Generator.
	 * @ingroup primes
	 * @ingroup randiter
	 *
	 * Generates prime of specified length.
	 * @internal
	 * It is given by <code>nextprime(2^_bits-p)</code> where <code>size(p) < _bits</code>.
	 * @todo
	 * One could use Integer::random_exact.
	 */
	//template<class Prime_Type = integer>
	class RandomPrimeIterator {
	private:

		uint64_t 	_bits;  //!< common lenght of all primes
		integer        _shift;  //!< @internal used to set proper bit size
		integer        _prime;  //!< the generated prime.
		Givaro::IntPrimeDom _IPD; //!< empty struct dealing with primality.

	public:
		/*! Constructor.
		 * @param bits size of primes (in bits). Default is 27 so it
		 * can fit in a <code>Linbox::Modular<double></code>.
		 * @param seed if \c 0 a seed will be generated, otherwise, the
		 * provided seed will be use.
		 */
		RandomPrimeIterator(uint64_t bits = 27, uint64_t seed = 0) :
			_bits(bits), _shift(integer(1)<<_bits)
		{
			linbox_check(bits >1);
			if (! seed)
				setSeed( BaseTimer::seed() );
			else
				setSeed( seed );

			integer::random(_prime,_bits-1);
			_prime = _shift - _prime;
			_IPD.nextprimein(_prime);
		}

		typedef integer Prime_Type ;

		/** @brief operator++()
		 *  creates a new random prime.
		 */
		inline RandomPrimeIterator &operator ++ ()
		{
			integer::random(_prime,_bits-1);
			_prime = _shift - _prime;
			_IPD.nextprimein(_prime);
			return *this;
		}

		/** @brief get the random prime.
		 *  returns the actual prime.
		 */
		const Prime_Type &operator *  () const
		{
			return _prime;
		}

		/** @brief get the random prime.
		 *  returns the actual prime.
		 *  @warning a new prime is not generated.
		 */
		const Prime_Type & randomPrime() const
		{
			return _prime;
		}

		/** @brief Sets the seed.
		 *  Set the random seed to be \p ul.
		 *  @param ul the new seed.
		 */
		void static setSeed(uint64_t ul)
		{
			integer::seeding(ul);
		}

		void setBits(uint64_t bits) {
			_bits = bits;
			_shift=(integer(1)<<_bits);

			linbox_check(bits >1);

			integer::random(_prime,_bits-1);
			_prime = _shift - _prime;
			_IPD.nextprimein(_prime);

		}

		template<class _ModField>
		void setBitsField()
		{
			// std::cout << _bits << std::endl;
			integer k = FieldTraits<_ModField >::maxModulus();
			// std::cout << k << std::endl;
			uint64_t bits = (uint64_t)(k.bitsize());
			if (!bits)
				throw("weird");
			--bits;
			// std::cout << bits << std::endl;
			if (bits < _bits)
				setBits(bits);
			// std::cout << _bits << std::endl;
		}

		template<class _ModField>
		bool setBitsDelayedField(size_t n)
		{
			integer k = FieldTraits<_ModField >::maxModulus();
			int bits = (int)(k.bitsize());
	//std::cout << "maxmodbits: " << bits << std::endl;
			if (!bits) throw("weird");
			--bits;
			bits -= (int)((log((double)n)/2./M_LN2));
	//std::cout << "delcorrect: " << bits << std::endl;
			if (bits < 0) return false;
			if (bits < (int64_t)_bits) setBits((uint64_t)bits);
			return true;
		}


	};


	/*! @brief Random Prime Iterator.
	 * @ingroup primes
	 * @ingroup randiter
	 *
	 * Generates prime of size smaller than a prescribed one.
	 * This class is closer to the LinBox::RandIterArchetype.
	 * @todo
	 * one could create the same one on a LinBox::PID_double ?
	 */
	class RandomPrimeIter {

		uint64_t 	_bits;  //!< max length for all primes
		integer         _seed;  //!< the generated prime.
		Givaro::GivRandom _generator;
		Givaro::IntPrimeDom _IPD;

	public:
		/*! Constructor.
		 * @param bits max size of primes (in bits). Default is 30 so it
		 * can fit in a <code>Linbox::Modular<double></code>.
		 * @param seed if \c 0 a seed will be generated, otherwise, the
		 * provided seed will be use.
		 */
		RandomPrimeIter(uint64_t bits = 30, uint64_t seed = 0) :
			_bits(bits), _generator(seed)
		{
			linbox_check(bits >1);
			if (! seed)
				seed = BaseTimer::seed() ;
			else
				_seed = seed ;

			integer::seeding(_seed);
		}

		/// destructor.
		~RandomPrimeIter() {}

		/// prime type
		typedef integer Prime_Type ;
		/// copy constructor.
		/// @param R random iterator to be copied.
		RandomPrimeIter (const RandomPrimeIter &R) :
			_bits(R._bits), _seed(R._seed), _generator(R._generator)
		{}

		typedef integer Element ;

		/// copy.
		/// @param R random iterator to be copied.
		RandomPrimeIter &operator=(const RandomPrimeIter &R)
		{
			if (this != &R) {
				_bits = R._bits;
				_seed = R._seed;
			}
			return *this;
		}

		/** @brief get a random prime of maximum size \c _bits .
		 * @param[out] a a prime number
		 */
		integer & random (integer & a) const
		{
			integer::random(a,_bits);
			_IPD.nextprimein(a);

			while (a.bitsize()>_bits)
				_IPD.prevprimein(a);

			return a;
		}

		integer  random () const
		{
			integer a ;
			return random(a);
		}

		integer & random_exact (integer & a) const
		{
			// integer::random_exact(a,_bits);
			integer::random(a,_bits-1); //!@todo uses random_exact when givaro is released.
			a = (integer(1)<<_bits) - a;

			_IPD.nextprimein(a);
			while (a.bitsize()>_bits)
				_IPD.prevprimein(a);

			return a;
		}

		integer random_exact () const
		{
			integer a ;
			return random_exact(a);
		}

		integer & random_between (integer & a, uint64_t _low_bits) const
		{
			linbox_check(_low_bits < _bits);
			// integer::random_exact(a,_bits);
			uint64_t ze_bits = (_low_bits + ((_bits - _low_bits)*_generator())) ;
			linbox_check (!(ze_bits<_low_bits) && !(ze_bits>_bits));
			integer::random(a,ze_bits-1); //!@todo uses random_between when givaro is released.
			a = (integer(1)<<ze_bits) - a;

			_IPD.nextprimein(a);
			while (a.bitsize()>_bits)
				_IPD.prevprimein(a);

			linbox_check(a.bitsize() >= _low_bits && a.bitsize() <= _bits) ;

			return a;
		}

		integer random_between ( uint64_t _low_bits) const
		{
			// std::cout << "random between " << _low_bits << " and " << _bits << std::endl;
			integer a ;
			random_between(a,_low_bits);
			// std::cout << a << std::endl;
			return a ;
		}

		void setBits (uint64_t  bits)
		{
			_bits = bits;
		}

	};

}

#endif //__LINBOX_random_prime_iterator_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
