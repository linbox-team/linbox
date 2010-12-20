/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2007,2010 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr>
 * Modified by <brice.boyer@imag.fr> (RandomPrimeIter)
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*! @file randiter/random-prime.h
 * @ingroup randiter
 * @brief Generates random positive prime \ref integers.
 */

#ifndef __LINBOX_random_prime_iterator_H
#define __LINBOX_random_prime_iterator_H
#include <linbox/integer.h>
#include <linbox/util/timer.h>
#include "linbox/util/debug.h"

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

		unsigned int 	_bits;  //!< common lenght of all primes
		integer        _shift;  //!< @internal used to set proper bit size
		integer        _prime;  //!< the generated prime.

	public:
		/*! Constructor.
		 * @param bits size of primes (in bits)
		 * @param seed if \c 0 a seed will be generated, otherwise, the
		 * provided seed will be use.
		 */
		RandomPrimeIterator(unsigned int bits = 30, unsigned long seed = 0) :
			_bits(bits), _shift(integer(1)<<_bits)
		{
			linbox_check(bits >1);
			if (! seed)
				RandomPrimeIterator::setSeed( BaseTimer::seed() );
			else
				RandomPrimeIterator::setSeed( seed );

			integer::random(_prime,_bits-1);
			_prime = _shift - _prime;
			nextprime( _prime, _prime);
		}

		typedef integer Prime_Type ;

		/** @brief operator++()
		 *  creates a new random prime.
		 */
		inline RandomPrimeIterator &operator ++ ()
		{
			integer::random(_prime,_bits-1);
			_prime = _shift - _prime;
			nextprime( _prime, _prime);
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
		void static setSeed(unsigned long ul)
		{
			integer::seeding(ul);
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

		unsigned int 	_bits;  //!< max length for all primes
		integer         _seed;  //!< the generated prime.

	public:
		/*! Constructor.
		 * @param bits max size of primes (in bits)
		 * @param seed if \c 0 a seed will be generated, otherwise, the
		 * provided seed will be use.
		 */
		RandomPrimeIter(unsigned int bits = 30, unsigned long seed = 0) :
			_bits(bits)
		{
			linbox_check(bits >1);
			if (! seed)
				_seed = BaseTimer::seed() ;
			else
				_seed = seed ;

			integer::seeding(_seed);
		}

		/// destructor.
		~RandomPrimeIter() {}

		/// copy constructor.
		/// @param R random iterator to be copied.
		RandomPrimeIter (const RandomPrimeIter &R) :
			_bits(R._bits), _seed(R._seed)
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

		/** @brief get the random prime.
		 * @param[out] a the new prime number
		 */
		integer & random (integer & a) const
		{
			integer::random_exact(a,_bits);
			nextprime( a, a);
			while (a.bitsize()>_bits)
				prevprime(a,a);

			return a;
		}

		void setBits (unsigned int  bits)
		{
			_bits = bits;
		}

	};

}

#endif //__LINBOX_random_prime_iterator_H

