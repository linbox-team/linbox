/* linbox/randiter/generic.h
 * 2004 june, bds and Dan Roche starting from:
 * Copyright (C) 1999-2001 William J Turner,
 *               2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Changed LargeModularRandIter to ModularRandIter, parameterized on the
 * element type. This change is for compatibility with the changes in
 * field/modular.h
 *
 * Renamed from large-modular.h to modular.h
 * ------------------------------------
 * 2002-05-14 William J. Turner <wjturner@acm.org>
 *
 * Seeded random number generator in constructor.  _seed was never used
 * before.
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file randiter/generic.h
 * @ingroup randiter
 * @brief Genric random iterator.
 */
//#error "deprecated and not tested"

#ifndef __LINBOX_generic_randiter_H
#define __LINBOX_generic_randiter_H

#include <iostream>
#include <vector>

#include "time.h"
#include "linbox/integer.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/util/commentator.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox-config.h"

namespace LinBox
{

	/** Random field base element generator.
	 * This is a generator of random field elements that can be used with
	 * any field.  It initializes elements using rand().  For prime fields
	 * with p < 2^32, a near-uniform distrubution can be expected.  For
	 * larger fields or non-prime fields, a near-uniform distribution on an
	 * unspecified subset of the elements can be expected.
	 */
	template <class Field>
	class GenericRandIter {
	public:

		typedef typename Field::Element Element;

		/** Constructor from field, sampling size, and seed.
		 * The random field element iterator works in the field F, is seeded
		 * by seed, and it returns any one element with probability no more
		 * than <code>1/min (size, F.characteristic(c))</code>.
		 * A sampling size of zero means to sample from the entire prime subfield.
		 * A seed of zero means to use some arbitrary seed for the generator which will vary from run to run.
		 * @param F LinBox field in which to do arithmetic
		 * @param size constant integer reference of sample size from which to
		 *             sample (default = modulus of field)
		 * @param seed constant integer reference from which to seed random number
		 *             generator (default = 0)
		 */
		GenericRandIter (const Field &F,
				 const integer &size = 0,
				 const integer &seed = 0) :
			_field (F), _size (size), _seed (seed)
		{
			if (_seed == 0) _seed = time (NULL);

			integer cardinality;

			F.cardinality (cardinality);

			if ((_size == 0) || (_size > cardinality))
				_size = cardinality;

			linbox_check(cardinality>0); // could be -1

			commentator().report (10, INTERNAL_DESCRIPTION)
			<< "Created random generator with size " << _size
			<< " and seed " << _seed << std::endl;

			// Seed random number generator
			srand (_seed);
		}

		GenericRandIter (const GenericRandIter<Field> &R) :
			_field (R._field), _size (R._size), _seed (R._seed)
		{}

		~GenericRandIter () {}

		GenericRandIter<Field> &operator=(const GenericRandIter<Field> &R)
		{
			if (this != &R) { // guard against self-assignment
				_size = R._size;
				_seed = R._seed;
			}

			return *this;
		}

		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		Element &random (Element &a) const
		{ return _field.init(a,rand()); }

		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		ElementAbstract &random (ElementAbstract &a) const
		{
			Element tmp;

			random (tmp);
			return (a = ElementEnvelope <Field> (tmp));
		}

	private:

		/// Field in which arithmetic is done
		Field _field;

		/// Sampling size
		integer _size;

		/// Seed
		long _seed;

	}; // class GenericRandIter
}
#endif //__LINBOX_generic_randiter_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

