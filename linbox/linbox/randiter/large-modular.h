/* -*- mode: c; style: linux -*- */

/* linbox/randiter/large-modular.h
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

#ifndef __RANDITER_LARGE_MODULAR_H
#define __RANDITER_LARGE_MODULAR_H

#include <iostream>
#include <vector>

#include "time.h"
#include "linbox/integer.h"
#include "linbox/field/large-modular.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"

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
	class LargeModularRandIter
	{
	    public:

		/// element type
		typedef integer element;

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
		LargeModularRandIter (const LargeModular &F, 
				      const integer &size = 0, 
				      const integer &seed = 0)
			: _F (F), _size (size), _seed (seed)
		{ 
			if (_seed == 0) _seed = time (NULL);    

			integer cardinality; F.cardinality (cardinality);
			if ( (_size == 0) 
			     || ( (cardinality != integer (-1)) && (_size > cardinality) ) )
				_size = cardinality;
		}

		/** Copy constructor.
		 * Constructs LargeModularRandIter object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * @param  R LargeModularRandIter object.
		 */
		LargeModularRandIter (const LargeModularRandIter &R) 
			: _F (R._F), _size (R._size), _seed (R._seed) {}

		/** Destructor.
		 * This destructs the random field element generator object.
		 */
		~LargeModularRandIter () {}
    
		/** Assignment operator.
		 * Assigns LargeModularRandIter object R to generator.
		 * @param  R LargeModularRandIter object.
		 */
		LargeModularRandIter &operator=(const LargeModularRandIter &R)
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
		element &random (element &a) 
		{
			// Create new random elements
			integer temp_integer;
			integer card;
			temp_integer = static_cast<integer>((double (rand ())/RAND_MAX)*double (_size));
			temp_integer %= _F.cardinality (card);
			if (temp_integer < 0) temp_integer += card;
			return (a = temp_integer);
		}
 
		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		Element_abstract &random (Element_abstract &a) 
		{
			integer tmp;

			random (tmp);
			return (a = Element_envelope <LargeModular> (tmp));
		}

	    private:

		/// Field in which arithmetic is done
		LargeModular _F;

		/// Sampling size
		integer _size;
    
		/// Seed
		integer _seed;

	}; // class LargeModularRandIter

} // namespace LinBox 

#endif // _LARGE_MODULAR_RANDITER_
