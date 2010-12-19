/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
 * Written by <brice.boyer@imag.fr>
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

/*! @file randiter/random-integer.h
 * @ingroup randiter
 * @ingroup integer
 * @brief Generates random \ref integers.
 */

#ifndef __LINBOX_random_integer_iterator_H
#define __LINBOX_random_integer_iterator_H

#include "linbox/integer.h"
#include "linbox/util/timer.h"
#include "linbox/util/debug.h"

namespace LinBox
{

	/*!@brief Random Prime Generator.
	 * @ingroup integers
	 * @ingroup randiter
	 *
	 * Generates integers of specified length.
	 * @tparam _Unsigned if \c true, then only non negative integers
	 * are generated, if \c false, their sign is random.
	 */
	template<bool _Unsigned=true>
	class RandomIntegerIterator {

		unsigned int    _bits;  //!< common lenght of all integers
		integer      _integer;  //!< the generated integer.

	public:
		/*! Constructor.
		 * @param bits size of integers (in bits)
		 * @param seed if \c 0 a seed will be generated, otherwise, the
		 * provided seed will be use.
		 */
		RandomIntegerIterator(unsigned int bits = 30, unsigned long seed = 0) :
			_bits(bits)
		{
			linbox_check(bits>1);
			if (! seed)
				RandomIntegerIterator::setSeed( BaseTimer::seed() );
			else
				RandomIntegerIterator::setSeed( seed );

			integer::random_exact<_Unsigned>(_integer,_bits);
		}

		typedef integer Integer_Type ;

		/** @brief operator++()
		 *  creates a new random integer.
		 */
		inline RandomIntegerIterator &operator ++ ()
		{
			integer::random_exact<_Unsigned>(_integer,_bits);
			return *this;
		}

		/** @brief get the random integer.
		 *  returns the actual integer.
		 */
		const Integer_Type &operator *  () const
		{
			return _integer;
		}

		/** @brief get the random integer.
		 *  returns the actual integer.
		 *  @warning a new integer is not generated.
		 */
		const Integer_Type & randomInteger() const
		{
			return _integer;
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

	/*! @brief Random Integer Iterator.
	 * @ingroup integers
	 * @ingroup randiter
	 *
	 * Generates integers of size smaller than a prescribed one.
	 * This class is closer to the LinBox::RandIterArchetype.
	 * @todo
	 * one could create the same one on a LinBox::PID_double ?
	 * @tparam _Unsigned if \c true, then only non negative integers
	 * are generated, if \c false, their sign is random.
	 */
	template<bool _Unsigned=true>
	class RandomIntegerIter {

		unsigned int 	_bits;  //!< max length for all integers
		integer         _seed;  //!< the generated integer.

	public:
		/*! Constructor.
		 * @param bits max size of integers (in bits)
		 * @param seed if \c 0 a seed will be generated, otherwise, the
		 * provided seed will be use.
		 */
		RandomIntegerIter(unsigned int bits = 30, unsigned long seed = 0) :
			_bits(bits)
		{
			linbox_check(bits>1);
			if (! seed)
				_seed = BaseTimer::seed() ;
			else
				_seed = seed ;

			integer::seeding(_seed);
		}

		/// destructor.
		~RandomIntegerIter() {}

		/// copy constructor.
		/// @param R random iterator to be copied.
		RandomIntegerIter (const RandomIntegerIter &R) :
			_bits(R._bits), _seed(R._seed)
		{}

		typedef integer Element ;

		/// copy.
		/// @param R random iterator to be copied.
		RandomIntegerIter &operator=(const RandomIntegerIter &R)
		{
			if (this != &R) {
				_bits = R._bits;
				_seed = R._seed;
			}
			return *this;
		}

		/** @brief get the random integer.
		 * @param[out] a the new integer number
		 */
		const integer & random (integer & a) const
		{
			integer::random_exact<_Unsigned>(a,_bits);

			return a;
		}

		void setBits (unsigned int  bits)
		{
			_bits = bits;
		}

		unsigned int getBits () const
		{
			return _bits ;
		}

	};

}

#endif //__LINBOX_random_integer_iterator_H

