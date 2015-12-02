/* Copyright (C) 2010 LinBox
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 * (adapted form random-prime.h)
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

namespace std {
    template <bool B>
    using bool_constant = integral_constant<bool, B>;
}

namespace LinBox
{
	/*!@brief Random Integer Iterator.
	 * @ingroup integers
	 * @ingroup randiter
	 *
	 * Generates integers of specified length.
	 * @tparam _Unsigned if \c true, then only non negative integers
	 * are generated, if \c false, their sign is random.
	 */
	template<bool _Unsigned=true, bool _Exact_Size=false>
	class RandomIntegerIterator {

		size_t    _bits;  //!< common lenght of all integers
		integer   _integer;  //!< the generated integer.

        typedef typename 
        std::bool_constant<_Exact_Size>::type _Exact_Size_t;
        
	public:
		typedef integer Integer_Type ;
		typedef integer Element ;

		/*! Constructor.
		 * @param bits size of integers (in bits)
		 * @param seed if \c 0 a seed will be generated, otherwise, the
		 * provided seed will be use.
		 */
		RandomIntegerIterator(size_t bits = 30, uint64_t seed = 0) :
			_bits(bits)
		{
			linbox_check(bits>1);
			if (! seed)
				setSeed( static_cast<uint64_t>(BaseTimer::seed()) );
			else
				setSeed( seed );

			this->operator++();
		}

		/// copy constructor.
		/// @param R random iterator to be copied.
		RandomIntegerIterator (const RandomIntegerIterator &R) :
			_bits(R._bits), _integer(R._integer)
		{}

		/// copy.
		/// @param R random iterator to be copied.
		RandomIntegerIterator &operator=(const RandomIntegerIterator &R)
		{
			if (this != &R) {
				_bits = R._bits;
				_integer = R._integer;
			}
			return *this;
		}


		/** @brief operator++()
		 *  creates a new random integer.
		 */
		inline RandomIntegerIterator &operator ++ ()
		{
			this->nextRandom(_Exact_Size_t(), _integer);
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

		Integer_Type & random (Integer_Type & a) const
		{
			return this->nextRandom(_Exact_Size_t(), a);
		}

		/** @brief Sets the seed.
		 *  Set the random seed to be \p ul.
		 *  @param ul the new seed.
		 */
		void static setSeed(uint64_t ul)
		{
			integer::seeding(ul);
		}

		void setBits (size_t  bits)
		{
			_bits = bits;
		}

		size_t getBits () const
		{
			return _bits ;
		}

        const Givaro::ZRing<Integer_Type> ring() const
        { 
            return Givaro::ZRing<Integer_Type>();
        }
        
                

    protected:
        inline Integer_Type& nextRandom(std::true_type,  Integer_Type & a) const {
			return integer::random_exact<_Unsigned>(a,_bits);
        }
        inline Integer_Type& nextRandom(std::false_type, Integer_Type & a) const {
			return integer::random_lessthan<_Unsigned>(a,_bits);
        }


	};

}

#endif //__LINBOX_random_integer_iterator_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
