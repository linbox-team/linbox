/* linbox/randiter/modular-balanced.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@acm.org>,
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
 * 2005-06-24 William J. Turner <wjturner@acm.org>
 *
 * Removed using declarations.
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

#ifndef __LINBOX_randiter_modular_balanced_H
#define __LINBOX_randiter_modular_balanced_H

#include <iostream>
#include <vector>

#include "time.h"
#include "linbox/integer.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/util/commentator.h"
//#include "linbox/randiter/mersenne-twister.h"
#include "linbox/linbox-config.h"

namespace LinBox
{

#if 0
	 template <class Element>
	 class Givaro::ModularBalanced<Element>;

	/** Random field base element generator.
	 * This encapsulated class is a generator of random field base elements for
	 * the encapsulating field.
	 * It is required to contain constructors from a field object and
	 * two integers.  The first integer being a cardinality of a set to
	 * draw the random elements from, and the second being a seed for the
	 * random number generator.
	 * It is also required to contain a copy constructor, a destructor, and
	 * an operator () which acts on a reference to a field base element.  In this
	 * operator (), the random element is placed into the input field base element
	 * and also returned as a reference.
	 */
	template <class Element>
	class Givaro::ModularBalancedRandIter {
	public:

		/** Constructor from field, sampling size, and seed.
		 * The random field element iterator works in the field F, is seeded
		 * by seed, and it returns any one element with probability no more
		 * than 1/min (size, F.cardinality (c)).
		 * A sampling size of zero means to sample from the entire field.
		 * A seed of zero means to use some arbitrary seed for the generator.
		 * Purely virtual.
		 * @param F LinBox field archetype object in which to do arithmetic
		 * @param size constant integer reference of sample size from which to
		 *             sample (default = modulus of field)
		 * @param seed constant integer reference from which to seed random number
		 *             generator (default = 0)
		 */
		Givaro::ModularBalancedRandIter (const Givaro::ModularBalanced<Element> &F,
					 const integer &size = 0,
					 const unit64_t seed = 0) :
			_field (F), _size (size), _seed (seed)
		{
			if (_seed == 0) _seed = time (NULL);

			integer cardinality;

			F.cardinality (cardinality);

			if ((_size == 0) || (_size > cardinality))
				_size = cardinality;

			commentator().report (10, INTERNAL_DESCRIPTION)
			<< "Created random generator with size " << _size
			<< " and seed " << _seed << std::endl;

			// Seed random number generator
			srand ((unsigned)_seed);
		}

		/** Copy constructor.
		 * Constructs Givaro::ModularBalancedRandIter object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * @param  R Givaro::ModularBalancedRandIter object.
		 */
		Givaro::ModularBalancedRandIter (const Givaro::ModularBalancedRandIter<Element> &R) :
			_field (R._field), _size (R._size), _seed (R._seed)
		{}

		/** Destructor.
		 * This destructs the random field element generator object.
		 */
		~Givaro::ModularBalancedRandIter () {}

		/** Assignment operator.
		 * Assigns Givaro::ModularBalancedRandIter object R to generator.
		 * @param  R Givaro::ModularBalancedRandIter object.
		 */
		Givaro::ModularBalancedRandIter<Element> &operator=(const Givaro::ModularBalancedRandIter<Element> &R)
		{
			if (this != &R) { // guard against self-assignment
				_size = R._size;
				_seed = R._seed;
				_field = R._field;
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
		{
			return _field.init(a,rand()); }

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
			return (a = ElementEnvelope <Givaro::ModularBalanced<Element> > (tmp));
		}

	private:
		/// Field in which arithmetic is done
		Givaro::ModularBalanced<Element> _field;

		/// Sampling size
		integer _size;

		/// Seed
		uint64_t _seed;

	}; // class Givaro::ModularBalancedRandIter
#endif
} // namespace LinBox

#endif //__LINBOX_randiter_modular_balanced_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
