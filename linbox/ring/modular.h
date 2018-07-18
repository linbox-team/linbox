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
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
