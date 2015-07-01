/* -*- mode: c; style: linux -*- */

/* linbox/src/randiter/archetype.h
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

#ifndef __RANDITER_ARCHETYPE_H
#define __RANDITER_ARCHETYPE_H

#include "linbox/field/archetype.h"
#include "linbox/field/abstract.h"
#include "linbox/element/abstract.h"
#include "linbox/randiter/abstract.h"

namespace LinBox
{
	class Element_archetype;

	/** Random field element generator archetype.
	 * Archetype for the random field element generator
	 * common object interface to \Ref{LinBox}.
	 *
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
	class RandIter_archetype
	{
	    public:
    
		/** @name Common Object Interface.
		 * These methods are required of all \Ref{LinBox} field element generators.
		 */
		//@{
    
		/// Element type
		typedef Element_archetype Element;
    
		/** Constructor from field, sampling size, and seed.
		 * The random field element iterator works in the field F, is seeded
		 * by seed, and it returns any one element with probability no more
		 * than 1/min (size, F.cardinality (c)).
		 * A sampling size of zero means to sample from the entire field.
		 * A seed of zero means to use some arbitrary seed for the generator.
		 * In this implementation, this means copying the field to
		 * which F._field_ptr points, the element to which F._elem_ptr points, 
		 * and the random element generator to which F._randIter_ptr points.
		 * @param F LinBox field archetype object in which to do arithmetic
		 * @param size constant integer reference of sample size from which to 
		 *             sample (default = 0)
		 * @param seed constant integer reference from which to seed random number 
		 *             generator (default = 0)
		 */
		RandIter_archetype (const Field_archetype &F, 
				    const Integer &size = 0, 
				    const Integer &seed = 0)
			{ _randIter_ptr = F._randIter_ptr->construct (*F._field_ptr, size, seed); }

		/** Copy constructor.
		 * Constructs RandIter_archetype object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * In this implementation, this means copying the random field element
		 * generator to which R._randIter_ptr points.
		 * @param  R RandIter_archetype object.
		 */
		RandIter_archetype (const RandIter_archetype &R) 
			{ _randIter_ptr = R._randIter_ptr->clone (); }

		/** Destructor.
		 * This destructs the random field element generator object.
		 * In this implementation, this destroys the generator by deleting 
		 * the random generator object to which _randIter_ptr points.
		 */
		~RandIter_archetype () 
			{ delete _randIter_ptr; }
    
		/** Assignment operator.
		 * Assigns RandIter_archetype object R to generator.
		 * In this implementation, this means copying the generator to
		 * which R._randIter_ptr points.
		 * @param  R RandIter_archetype object.
		 */
		RandIter_archetype &operator= (const RandIter_archetype &R)
		{
			if (this != &R) { // guard against self-assignment
				if (_randIter_ptr != 0) delete _randIter_ptr;
				if (R._randIter_ptr != 0)_randIter_ptr = R._randIter_ptr->clone ();
			}
			return *this;
		}
 
		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * @return reference to random field element
		 */
		Element &random (Element &a)
			{ _randIter_ptr->random (*a._elem_ptr); return a; }

		//@} Common Object Iterface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all 
		 * \Ref{LinBox Random field element generators}
		 * and are included only for this implementation of the archetype.
		 */
		//@{
    
		/** Constructor.
		 * Constructs field from pointer to \Ref{RandIter_abstract}.
		 * Not part of the interface.
		 * Creates new copies of random iterator generator object in dynamic memory.
		 * @param  randIter_ptr  pointer to \Ref{RandIter_abstract}
		 */
		RandIter_archetype (RandIter_abstract* randIter_ptr)
			: _randIter_ptr (randIter_ptr->clone ()) {}

		//@} Implementation-Specific Methods
    
	    private:

		/** Pointer to RandIter_abstract object.
		 * Not part of the interface.
		 * Included to allow for archetype use three.
		 */
		mutable RandIter_abstract *_randIter_ptr;
     
	}; // class RandIter_archetype
 
} // namespace LinBox

#endif // __RANDITER_ARCHETYPE_H
