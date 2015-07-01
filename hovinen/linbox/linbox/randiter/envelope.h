/* -*- mode: c; style: linux -*- */

/* linbox/src/randiter/envelope.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __RANDITER_ENVELOPE_H
#define __RANDITER_ENVELOPE_H

#include <iostream>
#include "linbox/field/envelope.h"
#include "linbox/element/envelope.h"
#include "linbox/randiter/abstract.h"

namespace LinBox 
{ 

	/** Random field base element generator.
	 * This encapsulated class is a generator of random field base elements for 
	 * the encapsulating field.
	 * It is required to contain constructors from a field object and
	 * two integers.  The first integer being a cardinality of a set to 
	 * draw the random elements from, and the second being a seed for the 
	 * random number generator.
	 * It is also required to contain a copy constructor, a destructor, and
	 * an operator() which acts on a reference to a field base element.  In this 
	 * operator(), the random element is placed into the input field base element 
	 * and also returned as a reference.
	 */
	template <class Field>
	class RandIter_envelope : public RandIter_abstract
	{
	    public:

		/// Element type
		typedef Element_envelope<Field> Element;

		/** Constructor from field, sampling size, and seed.
		 * The random field element iterator works in the field F, is seeded
		 * by seed, and it returns any one element with probability no more
		 * than 1/min (size, F.cardinality (c)).
		 * A sampling size of zero means to sample from the entire field.
		 * A seed of zero means to use some arbitrary seed for the generator.
		 * @param F LinBox field envelope object in which to do arithmetic
		 * @param size constant integer reference of sample size from which to 
		 *             sample (default = 0)
		 * @param seed constant integer reference from which to seed random number
		 *             generator (default = 0)
		 */
		RandIter_envelope (const Field_envelope<Field> &F, 
				   const Integer &size = 0, 
				   const Integer &seed = 0)
			: _randIter (F._field, size, seed) {}

		/** Constructor from random field element generator to be wrapped
		 * @param R random field element generator object to be wrapped
		 */
		RandIter_envelope (const typename Field::randIter &R) : _randIter (R) {}

		/** Copy constructor.
		 * Constructs RandIter_envelope object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * @param  R RandIter_envelope object.
		 */
		RandIter_envelope (const RandIter_envelope &R) : _randIter (R._randIter) {}

		/** Destructor.
		 * Required by abstract base class.
		 * This destructs the random field element generator object.
		 */
		~RandIter_envelope () {}
    
		/** Assignment operator.
		 * Assigns RandIter_envelope object R to generator.
		 * Required by abstract base class.
		 * @param  R RandIter_envelope object.
		 */
		RandIter_abstract &operator= (const RandIter_abstract &R)
		{
			if (this != &R) // guard against self-assignment
				_randIter = static_cast<const RandIter_envelope&> (R)._randIter;

			return *this;
		}
 
		/** Virtual constructor from field, sampling size, and seed.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * The random field element iterator works in the field F, is seeded
		 * by seed, and it returns any one element with probability no more
		 * than 1/min (size, F.cardinality (c)).
		 * A sampling size of zero means to sample from the entire field.
		 * A seed of zero means to use some arbitrary seed for the generator.
		 * Required by abstract base class.
		 * @param F LinBox field abstract object in which to do arithmetic
		 * @param size constant integer reference of sample size from which to 
		 *             sample (default = 0)
		 * @param seed constant integer reference from which to seed random number
		 *             generator (default = 0)
		 */
		RandIter_abstract *construct (const Field_abstract &F, 
					      const Integer &size = 0, 
					      const Integer &seed = 0) const
		{ 
			return new RandIter_envelope (static_cast<const Field_envelope<Field>&> (F)._field, size, seed);
		}

		/** Virtual copy constructor.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * Required by abstract base class.
		 * @return pointer to new RandIter_abstract object in dynamic memory.
		 */
		RandIter_abstract* clone (void) const
			{ return new RandIter_envelope (*this); }

		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		Element_abstract &random (Element_abstract &a)
			{ return _randIter.random (a); }

	    private:

		typename Field::randIter _randIter;

	}; // class RandIter_envelope

} // namespace LinBox

#endif // __RANDITER_ENVELOPE_H

