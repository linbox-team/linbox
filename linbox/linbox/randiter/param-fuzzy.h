/* -*- mode: c; style: linux -*- */

/* linbox/src/randiter/param-fuzzy.h
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

#ifndef __RANDITER_PARAM_FUZZY_H
#define __RANDITER_PARAM_FUZZY_H

#include <iostream>
#include <vector>

#include "linbox/integer.h"
#include "linbox/field/param-fuzzy.h"

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
	 * an operator () which acts on a reference to a field base element.  In this 
	 * operator (), the random element is placed into the input field base element 
	 * and also returned as a reference.
	 */
	class ParamFuzzyRandIter
	{
	    public:

		/// Element type
		typedef double Element;

		/** Constructor from field, sampling size, and seed.
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
		ParamFuzzyRandIter (const ParamFuzzy &F, 
				   const integer &size = 0, 
				   const integer &seed = 0)
			: _F (F), _size (size), _seed (seed)
		{ 
			if (_size == 0) F.cardinality (_size);
			if (_seed == 0) _seed = time (NULL);    
		}

		/** Copy constructor.
		 * Constructs ParamFuzzyRandIter object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * @param  R ParamFuzzyRandIter object.
		 */
		ParamFuzzyRandIter (const ParamFuzzyRandIter &R) 
			: _F (R._F), _size (R._size), _seed (R._seed) {}

		/** Destructor.
		 * This destructs the random field element generator object.
		 */
		~ParamFuzzyRandIter () {}
    
		/** Assignment operator.
		 * Assigns ParamFuzzyRandIter object R to generator.
		 * @param  R ParamFuzzyRandIter object.
		 */
		ParamFuzzyRandIter &operator=(const ParamFuzzyRandIter &R)
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
		Element &random (Element &a) 
		{
			// Create new random elements
			if (_size == 0)
				return (a = Element (rand ()));
			else
				return (a = Element (static_cast<long>((double (rand ())/RAND_MAX)*double (_size))));
		}

	    private:

		/// Field in which arithmetic is done
		ParamFuzzy _F;

		/// Sampling size
		Integer _size;
    
		/// Seed
		Integer _seed;

	}; // class ParamFuzzyRandIter : public ParamFuzzyRandIter

} // namespace LinBox 

#endif // __PARAM_FUZZY_H
