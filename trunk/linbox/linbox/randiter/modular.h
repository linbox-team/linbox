/* -*- mode: c; style: linux -*- */

/* linbox/randiter/modular.h
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
 *
 * See COPYING for license information.
 */

#ifndef __RANDITER_MODULAR_H
#define __RANDITER_MODULAR_H

#include <iostream>
#include <vector>

#include "time.h"
#include "linbox/integer.h"
#include "linbox/field/modular.h"
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
	template <class element>
	class ModularRandIter
	{
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
		ModularRandIter (const Modular<element> &F, 
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
		 * Constructs ModularRandIter object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * @param  R ModularRandIter object.
		 */
		ModularRandIter (const ModularRandIter<element> &R) 
			: _F (R._F), _size (R._size), _seed (R._seed) {}

		/** Destructor.
		 * This destructs the random field element generator object.
		 */
		~ModularRandIter () {}
    
		/** Assignment operator.
		 * Assigns ModularRandIter object R to generator.
		 * @param  R ModularRandIter object.
		 */
		ModularRandIter<element> &operator=(const ModularRandIter<element> &R)
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
			element temp_integer;
			integer card;
			temp_integer = static_cast<element>((double (rand ())/RAND_MAX)*double (_size));
			temp_integer %= (element) _F.cardinality (card);
			if (temp_integer < 0) temp_integer += (element) card;
			return (a = temp_integer);
		}
 
		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		ElementAbstract &random (ElementAbstract &a) 
		{
			element tmp;

			random (tmp);
			return (a = ElementEnvelope <Modular<element> > (tmp));
		}

	    private:

		/// Field in which arithmetic is done
		Modular<element> _F;

		/// Sampling size
		integer _size;
    
		/// Seed
		integer _seed;

	}; // class ModularRandIter

	template <class element>
	class Modular<element>::RandIter {
		ModularRandIter<element> _r;

	    public:
		RandIter (const Modular<element> &F, const integer &size = 0, const integer &seed = 0)
			: _r (F, size, seed) {}
		RandIter (const Modular<element>::RandIter &r)
			: _r (r._r) {}
		~RandIter () {}
		RandIter &operator= (const RandIter &r)
			{ _r = r._r; }
		element &random (element &a)
			{ return _r.random (a); }
		ElementAbstract &random (ElementAbstract &a) 
			{ return _r.random (a); }
	};

} // namespace LinBox 

#endif // _LARGE_MODULAR_RANDITER_
