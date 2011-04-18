/* linbox/randiter/envelope.h
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

#ifndef __LINBOX_randiter_envelope_H
#define __LINBOX_randiter_envelope_H

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
	class RandIterEnvelope : public RandIterAbstract
	{
	    public:

		/// element type
                //typedef ElementAbstract element;
		typedef ElementEnvelope<Field> Element;

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
		RandIterEnvelope (const FieldEnvelope<Field> &F, 
				   const integer &size = 0, 
				   const integer &seed = 0)
			: _randIter (F._field, size, seed) {}

		/** Constructor from random field element generator to be wrapped
		 * @param R random field element generator object to be wrapped
		 */
		RandIterEnvelope (const typename Field::RandIter &R) : _randIter (R) {}

		/** Copy constructor.
		 * Constructs RandIterEnvelope object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * @param  R RandIterEnvelope object.
		 */
		RandIterEnvelope (const RandIterEnvelope &R) : _randIter (R._randIter) {}

		/** Destructor.
		 * Required by abstract base class.
		 * This destructs the random field element generator object.
		 */
		~RandIterEnvelope () {}
    
		/** Assignment operator.
		 * Assigns RandIterEnvelope object R to generator.
		 * Required by abstract base class.
		 * @param  R RandIterEnvelope object.
		 */
		RandIterAbstract &operator= (const RandIterAbstract &R)
		{
			if (this != &R) // guard against self-assignment
				_randIter = static_cast<const RandIterEnvelope&> (R)._randIter;

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
		RandIterAbstract *construct (const FieldAbstract &F, 
					      const integer &size = 0, 
					      const integer &seed = 0) const
		{ 
			return new RandIterEnvelope (static_cast<const FieldEnvelope<Field>&> (F)._field, size, seed);
		}

		/** Virtual copy constructor.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * Required by abstract base class.
		 * @return pointer to new RandIterAbstract object in dynamic memory.
		 */
		RandIterAbstract* clone (void) const
			{ return new RandIterEnvelope (*this); }

		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		ElementAbstract &random (ElementAbstract &a) const
			//{ return  _randIter.random (a); }
			// GV Thu Apr 18 14:46:46 MEST 2002
			// modify by P.G. 2004-07-16
                     {
			     _randIter.random(static_cast<ElementEnvelope<Field>&> (a)._elem );
			     return  a; 
		     }



	    private:

		typename Field::RandIter _randIter;

	}; // class RandIterEnvelope

} // namespace LinBox

#endif // __LINBOX_randiter_envelope_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
