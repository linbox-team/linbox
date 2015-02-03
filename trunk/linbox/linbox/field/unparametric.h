
/* linbox/field/unparametric.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 *               2005 Clement Pernet
 *
 * Written by W. J. Turner <wjturner@acm.org>,
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

#ifndef __LINBOX_field_unparametric_H
#define __LINBOX_field_unparametric_H

#include <typeinfo>
#include <string>
#include <algorithm>
#include <givaro/unparametric.h>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/randiter/nonzero.h"


namespace LinBox
{

	template <class Ring>
	struct ClassifyRing;

	template <class K>
	struct ClassifyRing<Givaro::UnparametricRing<K> > {
		typedef RingCategories::GenericTag categoryTag;
	};

	template <>
	struct ClassifyRing<Givaro::UnparametricRing<integer> > {
		typedef RingCategories::IntegerTag categoryTag;
	};


	template<class Field>
	class FieldAXPY;

	/*! @ingroup integers
	 * @brief NO DOc
	 */
	template<>
	class FieldAXPY<Givaro::UnparametricRing<integer> >  {
	public:
		typedef Givaro::UnparametricRing<integer> Field;
		typedef integer Element;
		typedef Element Abnormal;

		/** Constructor.
		 * A faxpy object if constructed from a Field and a field element.
		 * Copies of this objects are stored in the faxpy object.
		 * @param F field F in which arithmetic is done
		 */
		FieldAXPY (const Field &F) :
			_field (F)
		{ _y = 0; }

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
		 * allow optimal multiplication, such as integer * int
		 */
		template<class Element1>
		inline Element&  mulacc (const Element &a, const Element1 &x)
		{
			return _y += (a * x);
		}

		template<class Element1>
		inline Element& accumulate (const Element1 &t)
		{
			return _y += t;
		}

		/** Retrieve y
		 *
		 * Performs the delayed modding out if necessary
		 */
		inline Element &get (Element &y) { y = _y; return y; }

		/** Assign method.
		 * Stores new field element for arithmetic.
		 * @return reference to self
		 * @param y constant reference to element a
		 */
		inline FieldAXPY &assign (const Element& y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			_y = 0;
		}
		
		const Field& field() const {
			return _field;
		}

	private:

		/// Field in which arithmetic is done
		Field _field;

		/// Field element for arithmetic
		Element _y;

	};

} // namespace LinBox

#include "linbox/randiter/unparametric.h"

#endif // __LINBOX_field_unparametric_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

