/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/element/abstract.h
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

	/** @memo Abstract element base class. A technicality.
	 * @doc The element class of \Ref{FieldAbstract}.
	 * This element has no knowledge of the field to which it belongs. 
	 * All operations and functions requiring knolwedge of the field,
	 * such as addition and other arithmetic operations, are supplied
	 * by the field and not the element class.
	 */
	class ElementAbstract 
	{
	    public:
    
		/** Virtual copy constructor.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * Purely virtual.
		 * @return pointer to new ElementAbstract object in dynamic memory.
		 */
		virtual ElementAbstract *clone (void) const = 0;

		/** Assignment operator.
		 * Purely virtual.
		 * @param  x constant reference to ElementAbstract object
		 * @return reference to self
		 */
		virtual ElementAbstract &operator= (const ElementAbstract &x) = 0;

		/** Destructor.
		 */
		virtual ~ElementAbstract (void) {}

	    protected:

		/** Default Constructor.
		 * Required by derived classes, but protected because this class should
		 * never be constructed by itself.
		 */
		ElementAbstract (void) {}

	}; // class ElementAbstract

} // namespace LinBox

#endif // __ELEMENT_ABSTRACT_H
