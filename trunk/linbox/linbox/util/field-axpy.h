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
	 * A fieldAXPY is an accumulator, allowing to add a number of field elements or products 
	 * in an unnormalized state, delaying modular reduction as long as possible.
	 * This class containse a value y and wraps the operations <code>y = y + a * x</code> and 
	 * <code>y = y + a</code>. 
	 *
	 * This default instance does no optimization, no delayed modular reduction.
	 * Through the use of template specialization, objects of this type can
	 * be used to speed up operations such as vector dot product operations. 
	 * In particular, for finite fields, dividing by the modulus and taking the remainder 
	 * is expensive. In many cases, this can be postponed until the end of the
	 * dot product operation, thus vastly improving performance.
	 *
	 * FieldAXPY<Fld> is an assignable type.  It is constructed from a field instance.
	 * The methods are mulacc(), accumulate(), operator=(), reset(), get(), and field().
	 * Of a const instance you can access get() and field().
	 * [Note: get() may renormalize the value, but it remains constant as a field element.]
	 *
	 * @param Field \ref LinBox @link Fields field@endlink
	 */
	template <class Field>
	class FieldAXPY {
	    public:

		/// Definition of element type
		typedef typename Field::Element Element;
		typedef Element Abnormal; // the type of the unnormalized values

		/** Constructor.
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
		FieldAXPY<Field> &operator = (const FieldAXPY<Field> &faxpy)
			{
				_field = faxpy.field();
				_y = faxpy._y;
				return *this;
			}

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

		/// reset value to zero.
		inline void reset() {
			field().init(_y,0);
		}

		inline const Field& field() const { return *_field; }
	    protected:

		/// Field in which arithmetic is done
		const Field *_field;

		/// Accumulator, unnormalized field element. 
		Abnormal _y;

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

