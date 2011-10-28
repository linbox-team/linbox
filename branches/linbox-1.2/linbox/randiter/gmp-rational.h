/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/randiter/gmp-rational.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file randiter/gmp-rational.h
 * @ingroup randiter
 * @brief Generate a random rational with prescribed size.
 */

#ifndef __LINBOX_randiter_gmp_rational_H
#define __LINBOX_randiter_gmp_rational_H

#include "linbox/field/gmp-rational.h"
#include "linbox/element/gmp-rational.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/linbox-config.h"

#include <sys/time.h>
#include <stdlib.h>

namespace LinBox
{

	class GMPRationalRandIter {
	public:

		/// Element type.
		typedef GMPRationalElement Element;

		/*! Rational number random iterator.
		 * @ingroup rationals
		 * @ingroup randiter
		 *
		 * @param F rational field
		 * @param size ??.
		 * @param seed seed for random generator. By default it is 0
		 * and a new one is created.
		 */
		GMPRationalRandIter (const GMPRationalField &F,
				     const integer &size = 0,
				     const integer &seed = 0) :
			_F (F), _size (size), _seed (seed)
		{
			if (seed == 0)
				_seed = time (NULL);
		}



		/// Copy constructor
		/// @param R Random iterator to copy.
		GMPRationalRandIter (const GMPRationalRandIter& R) :
			_F (R._F), _size (R._size), _seed (R._seed)
		{}


#ifdef __LINBOX_XMLENABLED
		GMPRationalRandIter(LinBox::Reader &R) :
			_F(R.Down(1))
		{
			R.Up(1);
			if(!R.expectTagName("randiter")) return;
			if(!R.expectAttributeNum("seed", _seed) || !R.expectAttributeNum("size", _size)) return;

			if(_seed == 0) _seed = time( NULL);

			return;

		}
#endif


		/// Destructor.
		~GMPRationalRandIter()
		{}

		/*! Copy.
		* @param R Random iterator to copy.
		* @return a reference to this copy.
		*/
		GMPRationalRandIter& operator=(const GMPRationalRandIter& R)
		{
			if (this != &R) { // guard against self-assignment
				_F = R._F;
				_seed = R._seed;
				_size = R._size;
			}
			return *this;
		}

		/*! get a random number.
		 * @param[out] a a random rational.
		 */
		Element &random (Element &a)  const
		{
			int value = 0;

			if (_size == 0) {
			unsigned int s;
				s = _seed;

				value = rand_r (&s);

				mpz_set_si (mpq_numref (a.rep), value);

				do {
					value = rand_r (&s);
				} while (value == 0);

				const_cast<integer&>(_seed) = s;
				mpz_set_si (mpq_denref (a.rep), value);
			}
			else {
				unsigned int s;
				int num, den;

				s = _seed;
				num = rand_r (&s);

				if (_size > 0) {
					int tmp = (int)_size;
					num %= tmp;
					den = 1L;
				}
				else {
					den = rand_r (&s);
				}

				const_cast<integer&>(_seed) = s;

				mpz_set_si (mpq_numref (a.rep), num);
				mpz_set_si (mpq_denref (a.rep), den);
			}

			mpq_canonicalize (a.rep);

			return a;
		}

		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @param[out] a random field element.
		 * @return reference to random field element \p a.
		 */
		ElementAbstract &random (ElementAbstract &a)  const
		{
			Element tmp;

			random (tmp);
			return (a = ElementEnvelope <GMPRationalField> (tmp));
		}

	private:

		GMPRationalField _F; //!< @internal rational field
		integer _size;       //!< @internal ??
		integer _seed;       //!< @internal seed

	}; // class GMPRationalRandIter

} // namespace LinBox

#endif // __LINBOX_randiter_gmp_random_H

