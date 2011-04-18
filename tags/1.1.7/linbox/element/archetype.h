/* linbox/element/archetype.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
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

#ifndef __LINBOX_element_archetype_H
#define __LINBOX_element_archetype_H

#include "linbox/element/abstract.h"

namespace LinBox
{

	// Forward declarations
	class FieldArchetype;
	class RandIterArchetype;

	/** @brief Field and Ring element interface specification and archetypical instance class.
\ingroup element

	 * Element classes must contain public default constructor, copy constructor,
	 * assignment operator, and destructor.  Note that primitive types
	 * such as int and double meet this specification.
	  
	 * In particular, the copy constructor is used to allow elements to 
	 * be passed by value to a function.
	  
	 * Note that the documentation below of the specific methods, describes
	 * them as implemented in the archetypic element class.

	 */
	class ElementArchetype
	{
	    public:

		/** @name Common Object Interface for LinBox Field elements.
		 * These methods are required of all \ref{LinBox} 
		 * {@link Fields field} elements.
		 */
		//@{

		/** Default constructor.
		 * This constructor is required to allow 
		 * {@link Fields field} elements to be primitive C++ types.
		 * Because constructor does not know what {@link Fields field} 
		 * the element belongs to, it cannot actually construct the element.
		 * In this implementation, the constructor it sets {\tt \_elem\_ptr}
		 * to the null pointer.  Initialization of the element is done through
		 * the field function init where the field is known.
		 */
		ElementArchetype (void) { _elem_ptr = 0; }

		/** Copy constructor.
		 * This constructor is required to allow 
		 * {@link Fields field} elements to be primitive C++ types, 
		 * and to allow field elements to be passed by value into 
		 * functions.
		 * Constructs {@link Fields field} element by copying the 
		 * {@link Fields field} element.
		 * In this implementation, this means copying the element to
		 * which {\tt a.\_elem\_ptr} points.
		 * @param  a field element.
		 */
		ElementArchetype (const ElementArchetype &a) 
		{ 
			if (a._elem_ptr != 0)
				_elem_ptr = a._elem_ptr->clone (); 
			else
				_elem_ptr = 0;
		}

		/** Destructor.
		 * In this implementation, this destroys element by deleting field 
		 * element to which {\tt \_elem\_ptr} points.
		 */
		~ElementArchetype () { if (_elem_ptr != 0) delete _elem_ptr; }

		/** Assignment operator.
		 * Assigns element a to element.  
		 * In this implementation, this is done 
		 * by copying field element to which {\tt \_elem\_ptr} points.
		 * @param  a field element.
		 */
		ElementArchetype &operator=(const ElementArchetype &a)
		{
			if (this != &a) { // guard against self-assignment
				if (_elem_ptr != 0) delete _elem_ptr;
				if (a._elem_ptr != 0) _elem_ptr = a._elem_ptr->clone ();
			}

			return *this;
		}

		//@} Common Object Interface

		/** @name Implementation-Specific Methods.
		 * These methods are not required of all LinBox field elements
		 * and are included only for this implementation of the archetype.
		 */
		//@{

		/** Constructor.
		 * Constructs field element from pointer to \ref{ElementAbstract}
		 * Not part of the interface.
		 * Creates new copy of element object in dynamic memory.
		 * @param  elem\_ptr  pointer to \ref{ElementAbstract}
		 */
		ElementArchetype (ElementAbstract *elem_ptr)
			: _elem_ptr (elem_ptr->clone ()) {}

		//@}
    
	    private:

		friend class RingArchetype;
		friend class FieldArchetype;
		friend class RandIterArchetype;
    
		/** @name Implementation-Specific Data.
		 * This data is not required of all LinBox field elements
		 * and is included only for this implementation of the archetype.
		 */
		//@{
    
		/** Pointer to parameterized field element.
		 * Not part of the common object interface for \ref{LinBox} field elements.
		 * Included to avoid code bloat.
		 */
		mutable ElementAbstract *_elem_ptr;
    
		//@} Non-Interface

	}; // class element

} // namespace LinBox

#endif // __LINBOX_element_archetype_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
