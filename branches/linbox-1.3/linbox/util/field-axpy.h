/* linbox/util/field-axpy.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_util_field_axpy_H
#define __LINBOX_util_field_axpy_H

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** FieldAXPY object.
	 *
	 * This class is used to wrap the operation <code>y = y + a * x</code>. It acts as an
	 * accumulator for \c y.
	 *
	 * Through the use of template specialization, objects of this type can
	 * be used to speed up vector dot product operations. In particular, for
	 * finite fields, dividing by the modulus and taking the remainder is
	 * expensive. In many cases, this can be postponed until the end of the
	 * dot product operation, thus vastly improving performance.
	 *
	 * This object is constructed from the field object \c F and a field
	 * element a which it stores and thus can use several times.  The use of
	 * an object instead of a static variable to store the element a makes
	 * this method thread-safe. (?? -bds)
	 *
	 * FieldAXPY<Fld> is an assignable type.  
	 * The methods are mulacc(), accumulate(), assign(), reset(), get(), and field().
	 * Of a const instance you can access get() and field().
	 * [Note: get() may renormalize the value, but it remains constant as a field element.
	 *
	 * @param Field \ref LinBox @link Fields field@endlink
	 */
	template <class Field>
	class FieldAXPY {
	    public:

		/// Definition of element type
		typedef typename Field::Element Element;

		/** Constructor.
		 * A faxpy object if constructed from a Field and a field element.
		 * Copies of this objects are stored in the faxpy object.
		 * @param F field F in which arithmetic is done
		 */
		FieldAXPY (const Field &F) :
		       	_field (&F)
		{ field().init (_y, 0); }

		/** Copy constructor.
		 * @param faxpy
		 */
		FieldAXPY (const FieldAXPY<Field> &faxpy) :
		       	_field (faxpy._field), _y (faxpy._y)
	       	{}

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
                { return field().axpyin (_y, a, x); }

            inline Element& accumulate (const Element &t)
                { return field().addin (_y, t); }

		/** Retrieve y
		 *
		 * Performs the delayed modding out if necessary
		 */
		inline Element &get (Element &y) const { return y = _y; }

		/** Assign method.
		 * Stores new field element for arithmetic.
		 * @return reference to self
		 * @param y constant reference to element a
		 */
		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			field().init(_y,0);
		}

		inline const Field& field() const { return *_field; }
	    protected:

		/// Field in which arithmetic is done
		/// Not sure why it must be mutable, but the compiler complains otherwise
		const Field *_field;

		/// Field element for arithmetic
		Element _y;

	}; // class FieldAXPY

} // namespace LinBox

#endif // __LINBOX_util_field_axpy_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

