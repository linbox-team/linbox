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
#include <givaro/modular-balanced.h>
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

	template <class Element, class Compute>
	struct ClassifyRing<Givaro::Modular<Element,Compute> const>
	{
		typedef RingCategories::ModularTag categoryTag;
	};

	template <class Element, class Compute>
	struct ClassifyRing<Givaro::Modular<Element,Compute>>
	{
		typedef RingCategories::ModularTag categoryTag;
	};

        template<class Element>
        struct ClassifyRing<Givaro::ModularBalanced<Element> > {
                typedef RingCategories::ModularTag categoryTag;
        };

        template<typename XXX>
        struct ClassifyRing<Givaro::GFqDom<XXX>> {
                typedef RingCategories::ModularTag categoryTag;
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


/*! Specialization of FieldAXPY and DotProducts for parameterized modular field */
#include "linbox/ring/modular/modular-int32.h"
#include "linbox/ring/modular/modular-int64.h"
#include "linbox/ring/modular/modular-short.h"
#include "linbox/ring/modular/modular-byte.h"
#include "linbox/ring/modular/modular-double.h"
#include "linbox/ring/modular/modular-float.h"
#include "linbox/ring/modular/modular-unsigned.h"

#endif // __LINBOX_field_modular_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

