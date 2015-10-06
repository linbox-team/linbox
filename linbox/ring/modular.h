/* linbox/ring/modular.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * LargeGivaro::Modular is now replace by a class Givaro::Modular parameterized on the element
 * type. So, the old LargeGivaro::Modular is equivalent to Givaro::Modular<integer>. All other
 * interface details are exactly the same.
 *
 * Renamed from large-modular.h to modular.h
 * ------------------------------------
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file ring/modular.h
 * @ingroup ring
 * @brief A Givaro::Modular ring is a representations of <code>Z/mZ</code>.
 * This file groups many implementations/specialisations of modular rings.
 *   - Givaro::Modular arithmetic is provided in the <code>ModularXXX<T></code> classes.
 *   - Specialisations for \ref FieldAXPY, \ref MVProductDomain, \ref DotProductDomain.
 *   - Random Iterators
 *   .
 *
 * @bug move Element& init(const Element&) to FFPACK. use using more..
 */

#ifndef __LINBOX_ring_modular_H
#define __LINBOX_ring_modular_H

#include <givaro/modular.h>
#include <givaro/gfq.h>
#include <iostream>
#include <climits>
#include <cmath>

#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/write-mm.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/linbox-config.h"
#include "linbox/field/field-traits.h"


// Namespace in which all LinBox code resides
namespace LinBox
{ /* Givaro::Modular Base */

	using Givaro::Caster;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<Givaro::Modular<Element> >
	{
		typedef RingCategories::ModularTag categoryTag;
	};

	template <class Element>
	struct ClassifyRing<Givaro::Modular<Element> const>
	{
		typedef RingCategories::ModularTag categoryTag;
	};

        template<typename XXX>
        struct ClassifyRing<Givaro::GFqDom<XXX>> {
                typedef RingCategories::ModularTag categoryTag;
        };




	/*! Specialization of FieldAXPY for parameterized modular field */

	template <class _Element>
	class FieldAXPY<Givaro::Modular<_Element> > {
	public:

		typedef _Element Element;
		typedef Element Abnormal;
		typedef Givaro::Modular<_Element> Field;

		FieldAXPY (const Field &F) :
			_field(&F), _y(F.zero)
		{}
		FieldAXPY (const FieldAXPY<Givaro::Modular<Element> > &faxpy) :
			_field (faxpy._field), _y (faxpy._y)
		{}

		FieldAXPY<Givaro::Modular <Element> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			return *this;
		}

		inline Element& mulacc (const Element &a, const Element &x)
		{
			return accumulate(a * x);
		}

		inline Element& accumulate (const Element &t)
		{
			return _y+=t;
		}

		inline Element &get (Element &y) { _y %= field().characteristic(); y = _y; return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			field().assign(_y, field().zero);
		}

		inline const Field &field() const { return *_field; }

	protected:

		const Field *_field;
		Element _y;
	};

	/*! Specialization of FieldAXPY for  modular double */

	template <>
	class FieldAXPY<Givaro::Modular<double,double> > {
	public:
        typedef Givaro::Modular<double,double> Field;
        typedef Field::Element Element;
		typedef Element Abnormal;

		FieldAXPY (const Field &F) :
			_field(F), _y(F.zero)
		{}
		FieldAXPY (const FieldAXPY<Givaro::Modular<Element> > &faxpy) :
			_field (faxpy._field), _y (faxpy._y)
		{}

		FieldAXPY<Field > &operator = (const FieldAXPY &faxpy)
		{
			const_cast<Field&>(_field) = faxpy._field;
			_y = faxpy._y;
			return *this;
		}

		inline Element& mulacc (const Element &a, const Element &x)
		{
			return accumulate(a * x);
		}

		inline Element& accumulate (const Element &t)
		{
			return _y+=t;
		}

		inline Element &get (Element &y) { return field().assign(y,field().reduce(_y)); }
        
        

		inline FieldAXPY &assign (const Element y)
		{
			field().assign(_y,y);
			return *this;
		}

		inline void reset() {
			field().assign(_y, field().zero);
		}

		inline const Field &field() const { return _field; }

	protected:

		const Field &_field;
		Element _y;
	};

	/*! Specialization of FieldAXPY for  modular float */

	template <>
	class FieldAXPY<Givaro::Modular<float,float> > {
	public:
        typedef Givaro::Modular<float,float> Field;
        typedef Field::Element Element;
		typedef Element Abnormal;

		FieldAXPY (const Field &F) :
			_field(F), _y(F.zero)
		{}
		FieldAXPY (const FieldAXPY<Givaro::Modular<Element> > &faxpy) :
			_field (faxpy._field), _y (faxpy._y)
		{}

		FieldAXPY<Field > &operator = (const FieldAXPY &faxpy)
		{
			const_cast<Field&>(_field) = faxpy._field;
			_y = faxpy._y;
			return *this;
		}

		inline Element& mulacc (const Element &a, const Element &x)
		{
			return accumulate(a * x);
		}

		inline Element& accumulate (const Element &t)
		{
			return _y+=t;
		}

		inline Element &get (Element &y) { return field().assign(y,field().reduce(_y)); }
        
        

		inline FieldAXPY &assign (const Element y)
		{
			field().assign(_y,y);
			return *this;
		}

		inline void reset() {
			field().assign(_y, field().zero);
		}

		inline const Field &field() const { return _field; }

	protected:

		const Field &_field;
		Element _y;
	};


/*
	template <>
	inline std::ostream& Givaro::ModularBase<Integer>::write (std::ostream &os) const
	{
		return os << "GMP integers mod " << _modulus;
	}
*/

	// template <>
	// inline integer& Givaro::Modular<integer>::init (integer& x, const double& y) const
	// {
	// 	integer tmp = (integer)y % _modulus;
	// 	if (tmp<0) tmp += _modulus;
	// 	return x = tmp;
	// }





} // namespace LinBox

#include "linbox/vector/vector-domain.h"

namespace LinBox {
	template<class Field>
	class MVProductDomain ;

} // LinBox


// #include "linbox/field/modular/modular-unsigned.h"
// #include "linbox/field/modular/modular-int32.h"
// // #ifdef __LINBOX_HAVE_INT64
// #include "linbox/field/modular/modular-int64.h"
// // #endif
// #include "linbox/field/modular/modular-short.h"
// #include "linbox/field/modular/modular-byte.h"
// #include "linbox/field/modular/modular-double.h"
// #include "linbox/field/modular/modular-float.h"

#endif // __LINBOX_field_modular_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

