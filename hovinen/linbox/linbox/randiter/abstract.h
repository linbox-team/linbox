/* -*- mode: c; style: linux -*- */

/* linbox/src/randiter/abstract.h
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

#ifndef __RANDITER_ABSTRACT_H
#define __RANDITER_ABSTRACT_H

#include <iostream>
#include <linbox/integer.h>

namespace LinBox 
{ 
	// forward declarations
	class Field_abstract;
	class Element_abstract;

	/** Random field element generator.
	 * This encapsulated class is a generator of random field elements for 
	 * the encapsulating field.
	 * It is required to contain constructors from a field object and
	 * two integers.  The first integer being a cardinality of a set to 
	 * draw the random elements from, and the second being a seed for the 
	 * random number generator.
	 * It is also required to contain a copy constructor, a destructor, and
	 * an operator () which acts on a reference to a field element.  In this 
	 * operator (), the random element is placed into the input field element 
	 * and also returned as a reference.
	 */
	class RandIter_abstract
	{
	    public:

		typedef Element_abstract Element;

		/** Virtual constructor from field, sampling size, and seed.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * The random field element iterator works in the field F, is seeded
		 * by seed, and it returns any one element with probability no more
		 * than 1/min (size, F.cardinality (c)).
		 * A sampling size of zero means to sample from the entire field.
		 * A seed of zero means to use some arbitrary seed for the generator.
		 * Purely virtual.
		 * @param F LinBox field archetype object in which to do arithmetic
		 * @param size constant integer reference of sample size from which to 
		 *             sample (default = 0)
		 * @param seed constant integer reference from which to seed random number
		 *             generator (default = 0)
		 */
		virtual RandIter_abstract *construct (const Field_abstract &F, 
						      const Integer &size = 0, 
						      const Integer &seed = 0) const = 0;

		/** Virtual copy constructor.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * Purely virtual.
		 * @return pointer to new RandIter_abstract object in dynamic memory.
		 */
		virtual RandIter_abstract *clone (void) const = 0;

		/** Assignment operator.
		 * Purely virtual.
		 * @param  x constant reference to RandIter_abstract object
		 * @return reference to self
		 */
		virtual RandIter_abstract &operator= (const RandIter_abstract &x) = 0;

		/** Destructor.
		 */
		virtual ~RandIter_abstract (void) {}

		/** Random field element creator.
		 * Purely virtual.
		 * @return reference to Element_abstract object
		 */
		virtual Element &random (Element &a) = 0;

	    protected:

		/** Default constructor
		 * Required by derived classes, but protected because this class should
		 * never be constructed by itself.
		 */
		RandIter_abstract (void) {}

	}; // class RandIter_abstract

} // namespace LinBox

#endif // __RANDITER_ABSTRACT_H
