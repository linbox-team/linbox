/* linbox/randiter/archetype.h
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

#ifndef __LINBOX_randiter_archetype_H
#define __LINBOX_randiter_archetype_H

#include "linbox/field/archetype.h"
#include "linbox/field/abstract.h"
#include "linbox/element/abstract.h"
#include "linbox/randiter/abstract.h"

namespace LinBox
{
	class ElementArchetype;

	/** \brief Random field element generator archetype.
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
	class RandIterArchetype
	{
	    public:
    
		/** @name Common Object Interface.
		 * These methods are required of all \Ref{LinBox} field element generators.
		 */
		//@{
    
		/// element type
		typedef ElementArchetype Element;
    
		/** Constructor from field, sampling size, and seed.
		 * The random field element iterator works in the
		 * field F, is seeded by seed, and it returns any one
		 * element with probability no more than 1/min (size,
		 * F.cardinality (c)).
		 *
		 * A sampling size of zero means to sample from the
		 * entire field.  A seed of zero means to use some
		 * arbitrary seed for the generator.  In this
		 * implementation, this means copying the field to
		 * which {\tt F.\_field\_ptr} points, the element to
		 * which {\tt F.\_elem\_ptr} points, and the random
		 * element generator to which {\tt F.\_randIter\_ptr} points.
		 *
		 * @param F LinBox field archetype object in which to do arithmetic
		 * @param size constant integer reference of sample size from which to 
		 *             sample (default = 0)
		 * @param seed constant integer reference from which to seed random number 
		 *             generator (default = 0)
		 */
		RandIterArchetype (const FieldArchetype &F, 
				    const integer &size = 0, 
				    const integer &seed = 0)
			{ _randIter_ptr = F._randIter_ptr->construct (*F._field_ptr, size, seed); }

		/** Copy constructor.
		 * Constructs RandIterArchetype object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * In this implementation, this means copying the random field element
		 * generator to which {\tt R.\_randIter\_ptr} points.
		 * @param  R RandIterArchetype object.
		 */
		RandIterArchetype (const RandIterArchetype &R) 
			{ _randIter_ptr = R._randIter_ptr->clone (); }


		/** Constructor.
		 * Constructs RandIterArchetype  from ANYTHING matching the interface
		 * using the enveloppe as a \Ref{FieldAbstract} and its
		 * encapsulated element and random element generator if needed.
		 * @param  field\_ptr pointer to field matching the interface
		 * @param  elem\_ptr  pointer to element matching the interface
		 * @param  randIter\_ptr  pointer to random matching the interface
		 */
		template<class Field_qcq>
			RandIterArchetype (Field_qcq *f, 
				    const integer &size = 0, 
				    const integer &seed = 0) 
		{ constructor (f, f, size, seed); }
	


		/** Destructor.
		 *
		 * This destructs the random field element generator
		 * object.  In this implementation, this destroys the
		 * generator by deleting the random generator object
		 * to which {\tt \_randIter\_ptr} points.
		 */
		~RandIterArchetype () 
			{ delete _randIter_ptr; }
    
		/** Assignment operator.
		 *
		 * Assigns RandIterArchetype object R to generator.
		 * In this implementation, this means copying the
		 * generator to which {\tt R.\_randIter\_ptr} points.
		 *
		 * @param  R RandIterArchetype object.
		 */
		RandIterArchetype &operator= (const RandIterArchetype &R)
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
		Element &random (Element &a) const
		{
			_randIter_ptr->random (*a._elem_ptr);
			return a;
		}

		//@} Common Object Iterface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all 
		 * \Ref{LinBox Random field element generators}
		 * and are included only for this implementation of the archetype.
		 */
		//@{
    
		/** Constructor.
		 * Constructs field from pointer to \Ref{RandIterAbstract}.
		 * Not part of the interface.
		 * Creates new copies of random iterator generator object in dynamic memory.
		 * @param  randIter_ptr  pointer to \Ref{RandIterAbstract}
		 */
		RandIterArchetype (RandIterAbstract* randIter_ptr)
			: _randIter_ptr (randIter_ptr->clone ()) {}

		//@} Implementation-Specific Methods
    
	    private:

		/** Pointer to RandIterAbstract object.
		 * Not part of the interface.
		 * Included to allow for archetype use three.
		 */
		mutable RandIterAbstract *_randIter_ptr;

		/** Template method for constructing archetype from a derived class of 
		 * FieldAbstract.
		 * This class is needed to help the constructor differentiate between 
		 * classes derived from FieldAbstract and classes that aren't.
		 * Should be called with the same argument to both parameters?
		 * @param	trait	pointer to FieldAbstract or class derived from it
		 * @param	field_ptr	pointer to class derived from FieldAbstract
		 */
		template<class Field_qcq>
		void constructor (FieldAbstract *trait, 
				  Field_qcq     *field_ptr, 
				  const integer &size = 0, 
				  const integer &seed = 0)
		{
			typename Field_qcq::RandIter FRI(*field_ptr);
			_randIter_ptr = FRI->construct (*field_ptr, size, seed);
		}
	 
		/** Template method for constructing archetype from a class not derived 
		 * from FieldAbstract.
		 * This class is needed to help the constructor differentiate between 
		 * classes derived from FieldAbstract and classes that aren't.
		 * Should be called with the same argument to both parameters?
		 * @param	trait	pointer to class not derived from FieldAbstract
		 * @param	field\_ptr	pointer to class not derived from FieldAbstract
		 */
		template<class Field_qcq>
		void constructor (void      *trait, 
				  Field_qcq *field_ptr, 
				  const integer &size = 0, 
				  const integer &seed = 0)
		{
			FieldEnvelope< Field_qcq > EnvF (*field_ptr);
			constructor (static_cast<FieldAbstract*> (&EnvF), &EnvF, size, seed) ;
		}

	}; // class RandIterArchetype
 
} // namespace LinBox

#endif // __LINBOX_randiter_archetype_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
