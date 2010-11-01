/* linbox/util/field-axpy.h
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

#ifndef __LINBOX_util_field_axpy_H
#define __LINBOX_util_field_axpy_H

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** FieldAXPY object.
	 *
	 * This class is used to wrap the operation y = y + a * x. It acts as an
	 * accumulator for y.
	 *
	 * Through the use of template specialization, objects of this type can
	 * be used to speed up vector dot product operations. In particular, for
	 * finite fields, dividing by the modulus and taking the remainder is
	 * expensive. In many cases, this can be postponed until the end of the
	 * dot product operation, thus vastly improving performance.
	 *
	 * This object is constructed from the field object F and a field
	 * element a which it stores and thus can use several times.  The use of
	 * an object instead of a static variable to store the element a makes
	 * this method thread-safe.
	 *
	 * @param Field \Ref{LinBox} {@link Fields field}
	 */
	template <class Field>
	class FieldAXPY 
	{
	    public:

		/// Definition of element type
		typedef typename Field::Element Element;

		/** Constructor.
		 * A faxpy object if constructed from a Field and a field element.
		 * Copies of this objects are stored in the faxpy object.
		 * @param F field F in which arithmetic is done
		 */
		FieldAXPY (const Field &F) : _F (F) { _F.init (_y, 0); }

		/** Copy constructor.
		 * @param faxpy
		 */
		FieldAXPY (const FieldAXPY<Field> &faxpy) : _F (faxpy._F), _y (faxpy._y) {}

		/** Assignment operator
		 * @param faxpy
		 */
		FieldAXPY<Field> &operator = (const FieldAXPY &faxpy) 
			{ _y = faxpy._y; return *this; }

		/** Add a*x to y
		 * y += a*x.
		 * @param a constant reference to element a
		 * @param x constant reference to element x
		 */
            inline Element& mulacc (const Element &a, const Element &x)
                { return _F.axpyin (_y, a, x); }
            
            inline Element& accumulate (const Element &t)
                { return _F.addin (_y, t); }
            
		/** Retrieve y
		 *
		 * Performs the delayed modding out if necessary
		 */
		inline Element &get (Element &y) { y = _y; return y; }

		/** Assign method.
		 * Stores new field element for arithmetic.
		 * @return reference to self
		 * @param y_init constant reference to element a
		 */
		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}
		
		inline void reset() {
			_F.init(_y,0);
		}

	    protected:

		/// Field in which arithmetic is done
		/// Not sure why it must be mutable, but the compiler complains otherwise
		Field _F;

		/// Field element for arithmetic
		Element _y;

	}; // class FieldAXPY

} // namespace LinBox

#endif // __LINBOX_util_field_axpy_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
