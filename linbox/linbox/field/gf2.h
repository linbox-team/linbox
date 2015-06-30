/* linbox/field/gf2.h
 * Copyright (C) 2003-2007 The LinBox group
 *
 * Authors : B. Hovinen, JG Dumas, C. Pernet
 * Moved into Givaro by : A. Breust
 *
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

#ifndef __LINBOX_field_gf2_H
#define __LINBOX_field_gf2_H

#include "linbox/field/field-traits.h"

#include <givaro/gf2.h>

namespace LinBox
{
	using GF2 = Givaro::GF2;

	template <class Ring>
	struct ClassifyRing;

	template<>
	struct ClassifyRing<GF2> {
		typedef RingCategories::ModularTag categoryTag;
	};
}

#endif

