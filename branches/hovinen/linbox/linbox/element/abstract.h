/* -*- mode: c; style: linux -*- */

/* linbox/src/element/abstract.h
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

#ifndef __ELEMENT_ABSTRACT_H
#define __ELEMENT_ABSTRACT_H

namespace LinBox 
{ 

	/** Abstract element base class.
	 * Element class of \Ref{Field_abstract}.
	 * This element has no knowledge of the field to which it belongs, 
	 * so all operations and functions requiring knolwedge of the field,
	 * such as addition and other arithmetic operations, must be supplied
	 * by the field and not the element.
	 */
	class Element_abstract 
	{
	    public:
    
		/** Virtual copy constructor.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * Purely virtual.
		 * @return pointer to new Element_abstract object in dynamic memory.
		 */
		virtual Element_abstract *clone (void) const = 0;

		/** Assignment operator.
		 * Purely virtual.
		 * @param  x constant reference to Element_abstract object
		 * @return reference to self
		 */
		virtual Element_abstract &operator= (const Element_abstract &x) = 0;

		/** Destructor.
		 */
		virtual ~Element_abstract (void) {}

	    protected:

		/** Default Constructor.
		 * Required by derived classes, but protected because this class should
		 * never be constructed by itself.
		 */
		Element_abstract (void) {}

	}; // class Element_abstract

} // namespace LinBox

#endif // __ELEMENT_ABSTRACT_H

